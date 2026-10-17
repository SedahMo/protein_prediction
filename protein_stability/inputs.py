"""Validate protein sequences and prepare explicitly mapped PDB coordinate files.

This module performs format and identity checks, never structure prediction or
stability scoring. PDB preparation keeps existing ATOM coordinates and numbering;
it does not repair missing atoms, build loops, add hydrogens, or minimize a model.
"""

from __future__ import annotations

import hashlib
import math
import os
from pathlib import Path
import re
import tempfile
from typing import Any


CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
RESIDUE_CODES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})
WATER_RESIDUES = frozenset({"HOH", "WAT", "DOD", "H2O"})
ALTLOC_POLICY = (
    "Keep blank/shared atoms and conformer A. Discard other conformers only when "
    "the same atom has conformer A. Reject duplicate atoms, blank/A conflicts, "
    "mixed residue identities, and atoms available only in another conformer."
)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read wrapped FASTA with unique first-token IDs and canonical amino acids.

    Lowercase sequences are normalized. Blank lines and outer line whitespace
    are accepted; internal whitespace, gaps, stop symbols, ambiguous residues,
    digits, comments, empty records, and duplicate identifiers are rejected.
    """
    records: dict[str, str] = {}
    current_id: str | None = None
    sequence_parts: list[str] = []

    def finish_record() -> None:
        if current_id is not None:
            if not sequence_parts:
                raise ValueError(f"FASTA record {current_id!r} has an empty sequence")
            records[current_id] = "".join(sequence_parts)

    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA line {line_number}: empty record header")
                current_id = header.split()[0]
                if current_id in records:
                    raise ValueError(f"FASTA line {line_number}: duplicate ID {current_id!r}")
                sequence_parts = []
                continue
            if current_id is None:
                raise ValueError(f"FASTA line {line_number}: sequence precedes the first header")
            sequence = line.upper()
            invalid = sorted(set(sequence) - CANONICAL_AA)
            if invalid:
                raise ValueError(
                    f"FASTA line {line_number}, record {current_id!r}: "
                    f"noncanonical sequence character(s) {', '.join(repr(c) for c in invalid)}"
                )
            sequence_parts.append(sequence)
    finish_record()
    if not records:
        raise ValueError("FASTA contains no protein records")
    return records


def inspect_fasta(path: str | Path) -> dict[str, Any]:
    """Return sequence identity/QC observations without a stability estimate."""
    sequences = read_fasta(path)
    records = []
    for identifier, sequence in sequences.items():
        warnings = []
        if sequence.count("C") % 2:
            warnings.append("Odd cysteine count; verify the intended disulfide design.")
        if re.search(r"(.)\1{5,}", sequence):
            warnings.append("Contains a homopolymer run of at least six residues.")
        if re.search(r"[GS]{12,}", sequence):
            warnings.append("Contains at least 12 consecutive glycine/serine residues; annotate the linker.")
        records.append({
            "id": identifier,
            "length": len(sequence),
            "sha256": _sha256(sequence.encode("ascii")),
            "cysteine_count": sequence.count("C"),
            "warnings": warnings,
        })
    return {
        "source": str(Path(path)),
        "source_sha256": _sha256(Path(path).read_bytes()),
        "record_count": len(records),
        "records": records,
        "interpretation": "Sequence QC only; no thermodynamic stability score is calculated.",
    }


def _validate_chain_id(chain: str) -> None:
    if not isinstance(chain, str) or not re.fullmatch(r"[A-Za-z0-9]", chain):
        raise ValueError("FoldX preparation requires single-character ASCII alphanumeric chain IDs")


def _parse_pdb(
    path: str | Path,
    selected_chains: set[str] | None,
    *,
    allow_insertions: bool,
) -> tuple[dict[str, Any], dict[str, list[str]]]:
    source = Path(path)
    raw_bytes = source.read_bytes()
    try:
        lines = raw_bytes.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise ValueError("PDB coordinate input must be ASCII fixed-width PDB, not mmCIF") from exc
    if any(line.startswith(("data_", "_atom_site.")) for line in lines[:100]):
        raise ValueError("mmCIF is unsupported; provide an explicitly converted fixed-width PDB")

    warnings: list[str] = []
    removed = {"water_atoms": 0, "ligand_atoms": 0,
               "unselected_protein_atoms": 0, "alternate_location_atoms": 0}
    # Each residue contains all conformers until a deterministic choice is made.
    residues: dict[tuple[str, int, str], dict[str, Any]] = {}
    chain_residues: dict[str, list[tuple[str, int, str]]] = {}
    last_residue: dict[str, tuple[str, int, str]] = {}
    last_chain: str | None = None
    last_atom_chain: str | None = None
    closed_chains: set[str] = set()
    terminated_chains: set[str] = set()
    model_count = 0
    inside_model = False
    saw_coordinate = False
    ended = False

    for line_number, line in enumerate(lines, 1):
        record = line[:6].strip()
        if record == "MODEL":
            model_count += 1
            if model_count > 1 or inside_model or saw_coordinate or ended:
                raise ValueError("PDB must contain exactly one coordinate model; split model ensembles first")
            inside_model = True
            continue
        if record == "ENDMDL":
            if not inside_model:
                raise ValueError(f"PDB line {line_number}: ENDMDL without an open MODEL")
            inside_model = False
            continue
        if record == "END":
            ended = True
            continue
        if record == "TER":
            # A short TER has no chain column: associate it with the last ATOM,
            # including an unselected chain, rather than the last selected one.
            ter_chain = line[21:22] if len(line) > 21 and line[21:22] != " " else last_atom_chain
            if ter_chain is not None:
                terminated_chains.add(ter_chain)
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        if ended or (model_count and not inside_model):
            raise ValueError(f"PDB line {line_number}: coordinate outside the single MODEL/END boundary")
        saw_coordinate = True
        if len(line) < 27:
            raise ValueError(f"PDB line {line_number}: truncated fixed-width coordinate record")
        resname = line[17:20].strip()
        if record == "HETATM":
            removed["water_atoms" if resname in WATER_RESIDUES else "ligand_atoms"] += 1
            continue
        chain = line[21:22]
        last_atom_chain = chain
        if selected_chains is not None and chain not in selected_chains:
            removed["unselected_protein_atoms"] += 1
            continue
        if chain in terminated_chains:
            raise ValueError(
                f"PDB chain {chain!r} reappears after TER; disconnected chain segments "
                "cannot be merged for FoldX execution"
            )
        if len(line) < 54:
            raise ValueError(f"PDB line {line_number}: ATOM coordinates require fixed-width columns 31–54")
        if resname not in RESIDUE_CODES:
            raise ValueError(
                f"PDB line {line_number}, chain {chain!r}: noncanonical residue {resname!r}; "
                "resolve chemical identity explicitly before FoldX"
            )
        try:
            number = int(line[22:26])
            serial = int(line[6:11])
            coordinates = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        except ValueError as exc:
            raise ValueError(f"PDB line {line_number}: invalid fixed-width atom serial, residue number, or xyz") from exc
        if not all(math.isfinite(value) for value in coordinates):
            raise ValueError(f"PDB line {line_number}: xyz coordinates must be finite")
        insertion = line[26:27]
        if insertion != " " and not allow_insertions:
            raise ValueError(
                f"PDB chain {chain!r}, residue {number}{insertion}: insertion codes are unsupported "
                "for FoldX execution; provide a deliberately remapped structure and FASTA map"
            )
        atom_name = line[12:16].strip()
        if not atom_name:
            raise ValueError(f"PDB line {line_number}: empty atom name")
        altloc = line[16:17]
        key = (chain, number, insertion)
        if last_chain != chain:
            if chain in closed_chains:
                raise ValueError(f"PDB chain {chain!r} reappears after another selected chain; resolve chain segments")
            if last_chain is not None:
                closed_chains.add(last_chain)
            last_chain = chain
        if key != last_residue.get(chain):
            if key in residues:
                raise ValueError(f"PDB contains a noncontiguous duplicate residue identifier: {chain}:{number}{insertion.strip()}")
            residues[key] = {"resname": resname, "atoms": {}, "serial": serial}
            chain_residues.setdefault(chain, []).append(key)
            last_residue[chain] = key
        residue = residues[key]
        if residue["resname"] != resname:
            raise ValueError(f"PDB residue {chain}:{number}{insertion.strip()} has ambiguous residue identities")
        conformers = residue["atoms"].setdefault(atom_name, {})
        if altloc in conformers:
            raise ValueError(
                f"PDB residue {chain}:{number}{insertion.strip()} has duplicate atom "
                f"{atom_name!r} with altloc {altloc!r}"
            )
        conformers[altloc] = (line, serial)

    if inside_model:
        raise ValueError("PDB has an unclosed MODEL; ENDMDL is required")
    if selected_chains is not None:
        missing = selected_chains - set(chain_residues)
        if missing:
            raise ValueError(f"Selected PDB chain(s) have no protein ATOM records: {', '.join(sorted(missing))}")
    if not residues:
        raise ValueError("PDB contains no selected canonical protein ATOM residues")

    chains: dict[str, Any] = {}
    residue_map: dict[str, str] = {}
    output_lines: dict[str, list[str]] = {}
    insertion_residues = []
    selected_atom_count = 0
    for chain, keys in chain_residues.items():
        sequence = []
        numbers = []
        chain_output = []
        residue_details = []
        previous_key: tuple[str, int, str] | None = None
        last_serial = 0
        for key in keys:
            _, number, insertion = key
            residue = residues[key]
            chosen_atoms = {}
            for atom_name, conformers in residue["atoms"].items():
                if " " in conformers:
                    if len(conformers) != 1:
                        raise ValueError(
                            f"PDB residue {chain}:{number}{insertion.strip()}, atom {atom_name}: "
                            "blank altloc conflicts with an alternate atom"
                        )
                    chosen_atoms[atom_name] = conformers[" "]
                elif "A" in conformers:
                    chosen_atoms[atom_name] = conformers["A"]
                    removed["alternate_location_atoms"] += len(conformers) - 1
                else:
                    raise ValueError(
                        f"PDB residue {chain}:{number}{insertion.strip()}, atom {atom_name}: "
                        "only non-A alternate locations are available"
                    )
            missing_atoms = BACKBONE_ATOMS - set(chosen_atoms)
            if missing_atoms:
                raise ValueError(
                    f"PDB residue {chain}:{number}{insertion.strip()} is missing backbone atom(s): "
                    f"{', '.join(sorted(missing_atoms))}; provide a completed structure"
                )
            # Keep original atom order, coordinates, chain IDs, residue numbering,
            # and atom serials. Clear the accepted A marker after conformer selection.
            for line, serial in chosen_atoms.values():
                chain_output.append((line[:16] + " " + line[17:]).rstrip())
                last_serial = serial
            selected_atom_count += len(chosen_atoms)
            amino_acid = RESIDUE_CODES[residue["resname"]]
            sequence.append(amino_acid)
            numbers.append(number)
            residue_id = f"{chain}:{number}{insertion.strip()}"
            residue_map[residue_id] = amino_acid
            residue_details.append({"number": number, "insertion_code": insertion.strip(),
                                    "amino_acid": amino_acid, "atom_names": list(chosen_atoms)})
            if insertion.strip():
                insertion_residues.append(residue_id)
            if previous_key is not None and number != previous_key[1] + 1:
                warnings.append(
                    f"Residue numbering discontinuity in chain {chain!r}: "
                    f"{previous_key[1]}{previous_key[2].strip()} -> {number}{insertion.strip()}; "
                    "no residues were inserted or renumbered."
                )
            previous_key = key
        final_key = keys[-1]
        final_resname = residues[final_key]["resname"]
        ter_serial = last_serial + 1
        if ter_serial > 99999:
            chain_output.append("TER")
        else:
            chain_output.append(
                f"TER   {ter_serial:5d}      {final_resname:>3} {chain}{final_key[1]:4d}{final_key[2]}"
            )
        output_lines[chain] = chain_output
        chains[chain] = {"sequence": "".join(sequence), "length": len(sequence),
                         "residue_numbers": numbers, "residues": residue_details}
    if removed["water_atoms"]:
        warnings.append(f"Excluded {removed['water_atoms']} water HETATM records.")
    if removed["ligand_atoms"]:
        warnings.append(
            f"Excluded {removed['ligand_atoms']} ligand/nonprotein HETATM records; "
            "assess whether any are required cofactors or covalent modifications."
        )
    if removed["unselected_protein_atoms"]:
        warnings.append(f"Excluded {removed['unselected_protein_atoms']} ATOM records from unselected chains.")
    if removed["alternate_location_atoms"]:
        warnings.append(
            f"Discarded {removed['alternate_location_atoms']} alternate-location atoms after choosing conformer A."
        )
    if insertion_residues:
        warnings.append("Insertion-coded residues are present; prepare_pdb will reject these for FoldX execution.")
    warnings.append(
        "Backbone presence and sequence identity are checked; side-chain completeness, geometry, "
        "protonation, disulfide connectivity, and model confidence require separate review."
    )
    return {
        "source": str(source), "source_sha256": _sha256(raw_bytes),
        "model_count": 1, "selected_chains": list(chains), "chains": chains,
        "residue_map": residue_map, "insertion_coded_residues": insertion_residues,
        "selected_atom_count": selected_atom_count, "removed": removed,
        "altloc_policy": ALTLOC_POLICY, "warnings": warnings,
    }, output_lines


def inspect_pdb(path: str | Path) -> dict[str, Any]:
    """Inspect all protein ATOM chains in one fixed-width PDB model.

    Insertion codes are retained for inspection. Preparation for FoldX rejects
    them. The same coordinate, residue, backbone, and altloc checks apply here.
    """
    audit, _ = _parse_pdb(path, None, allow_insertions=True)
    return audit


def prepare_pdb(
    source: str | Path,
    destination: str | Path,
    chain_map: dict[str, str],
    fasta: dict[str, str],
) -> dict[str, Any]:
    """Write selected protein atoms after exact chain-to-FASTA identity checks.

    ``chain_map`` maps each selected PDB chain ID to an existing FASTA record ID.
    Residue and atom numbering are retained. The output contains only selected
    ATOM records, one TER per chain, and END; accepted A altloc markers are cleared.
    The source file cannot be overwritten. No FoldX binary is needed.
    """
    if not isinstance(chain_map, dict) or not chain_map:
        raise ValueError("An explicit, nonempty PDB chain-to-FASTA-ID map is required")
    for chain, fasta_id in chain_map.items():
        _validate_chain_id(chain)
        if not isinstance(fasta_id, str) or fasta_id not in fasta:
            raise ValueError(f"PDB chain {chain!r} maps to missing FASTA ID {fasta_id!r}")
        sequence = fasta[fasta_id]
        if not isinstance(sequence, str) or not sequence or set(sequence) - CANONICAL_AA:
            raise ValueError(f"FASTA record {fasta_id!r} must be a nonempty uppercase canonical amino-acid sequence")
    source_path, destination_path = Path(source), Path(destination)
    if source_path.resolve() == destination_path.resolve():
        raise ValueError("PDB preparation destination must differ from the source file")
    audit, output_lines = _parse_pdb(source_path, set(chain_map), allow_insertions=False)
    for chain, data in audit["chains"].items():
        fasta_id = chain_map[chain]
        expected = fasta[fasta_id]
        observed = data["sequence"]
        if observed != expected:
            common_length = min(len(expected), len(observed))
            first_difference = next((i + 1 for i in range(common_length)
                                     if expected[i] != observed[i]), common_length + 1)
            raise ValueError(
                f"PDB chain {chain!r} sequence does not exactly match FASTA {fasta_id!r}: "
                f"ATOM length {len(observed)}, FASTA length {len(expected)}, "
                f"first difference at sequence position {first_difference}"
            )
        data["fasta_id"] = fasta_id
        data["sequence_sha256"] = _sha256(expected.encode("ascii"))
    output = "\n".join(line for chain in audit["selected_chains"]
                       for line in output_lines[chain]) + "\nEND\n"
    output_bytes = output.encode("ascii")
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile("wb", dir=destination_path.parent, delete=False) as handle:
            temporary_name = handle.name
            handle.write(output_bytes)
        os.replace(temporary_name, destination_path)
    finally:
        if temporary_name is not None and os.path.exists(temporary_name):
            os.unlink(temporary_name)
    audit["destination"] = str(destination_path)
    audit["output_sha256"] = _sha256(output_bytes)
    audit["chain_map"] = dict(chain_map)
    return audit
