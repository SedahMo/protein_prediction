"""Plan and execute licensed FoldX point-mutation studies on one backbone.

Nothing here predicts structures or bundles FoldX. Planning runs without FoldX.
Execution copies user-provided licensed resources into isolated working folders.
A declared version is recorded as such; the executable hash identifies the binary.
"""
from __future__ import annotations

import copy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
from typing import Any

from . import __version__
from .inputs import CANONICAL_AA, inspect_pdb, prepare_pdb, read_fasta
from .results import parse_diff


class FoldXError(ValueError):
    """A configuration or external FoldX execution failed validation."""


COMPARISON_REASONS = {
    "linker": "Linker sequence/length changes require independently validated conformational ensembles and matched experimental assays; BuildModel point substitutions cannot compare linker topologies.",
    "domain_order": "VH/VL order or connection changes alter topology; prepare validated models and compare geometry/developability plus matched experiments rather than rank BuildModel energies.",
    "cross_protein": "Different proteins/backbones do not share the mutant-specific WT reference required for comparable FoldX mutation ΔΔG; use matched within-backbone mutation studies and experimental stability assays.",
}
SOURCES = [
    "https://foldxsuite.crg.eu/command/RepairPDB",
    "https://foldxsuite.crg.eu/command/BuildModel",
    "https://foldxsuite.crg.eu/parameter/mutant-file",
    "https://foldxsuite.crg.eu/parameter/ionStrength",
]


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, data: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _object(value: Any, label: str, allowed: set[str]) -> dict:
    if not isinstance(value, dict):
        raise FoldXError(f"{label} must be an object")
    unknown = set(value) - allowed
    if unknown:
        raise FoldXError(f"{label}: unknown fields {sorted(unknown)}")
    return value


def _slug(value: Any, label: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", value) is None:
        raise FoldXError(f"{label} must be a 1–80 character identifier using letters, digits, _, - or .")
    return value


def _integer(value: Any, label: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise FoldXError(f"{label} must be an integer from {minimum} to {maximum}")
    return value


def _number(value: Any, label: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not minimum <= value <= maximum:
        raise FoldXError(f"{label} must be a finite number from {minimum} to {maximum}")
    return float(value)


def _file(value: Any, label: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise FoldXError(f"{label} requires a file path")
    path = Path(value).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FoldXError(f"{label} is missing, empty or not a regular file: {path}")
    return path


def validate_config(config: dict[str, Any]) -> dict[str, Any]:
    """Validate schema/types and mutation definitions without creating files.

    Residue numbers are original PDB numbers, not FASTA/IMGT numbering. This
    wrapper intentionally supports only positive PDB positions (1 through 9999).
    """
    _object(config, "config", {"schema_version", "study_id", "comparison_type", "synthetic", "reference", "variants", "foldx"})
    if type(config.get("schema_version")) is not int or config["schema_version"] != 1:
        raise FoldXError("schema_version must be integer 1")
    result = copy.deepcopy(config)
    result["study_id"] = _slug(config.get("study_id"), "study_id")
    comparison = config.get("comparison_type")
    if not isinstance(comparison, str):
        raise FoldXError("comparison_type must be a string")
    if comparison in COMPARISON_REASONS:
        raise FoldXError(COMPARISON_REASONS[comparison])
    if comparison != "point_mutation":
        raise FoldXError("comparison_type must be point_mutation, linker, domain_order or cross_protein")
    if not isinstance(config.get("synthetic", False), bool):
        raise FoldXError("synthetic must be a boolean")
    result["synthetic"] = config.get("synthetic", False)
    reference = _object(config.get("reference"), "reference", {"fasta", "pdb", "chain_map"})
    for name in ("fasta", "pdb"):
        result["reference"][name] = str(_file(reference.get(name), f"reference.{name}"))
    mapping = reference.get("chain_map")
    if not isinstance(mapping, dict) or not mapping:
        raise FoldXError("reference.chain_map must explicitly map selected PDB chains to FASTA IDs")
    for chain, fasta_id in mapping.items():
        if not isinstance(chain, str) or re.fullmatch(r"[A-Za-z0-9]", chain) is None:
            raise FoldXError("chain_map keys must be single alphanumeric PDB chain IDs")
        if not isinstance(fasta_id, str) or not fasta_id.strip():
            raise FoldXError("chain_map values must be nonempty FASTA record IDs")
    variants = config.get("variants")
    if not isinstance(variants, list) or not variants:
        raise FoldXError("variants must be a nonempty list of point-substitution designs")
    ids: set[str] = set()
    designs: set[tuple] = set()
    for index, variant in enumerate(variants):
        _object(variant, f"variants[{index}]", {"id", "mutations"})
        variant_id = _slug(variant.get("id"), f"variants[{index}].id")
        if variant_id in ids:
            raise FoldXError(f"Duplicate variant ID {variant_id}")
        ids.add(variant_id)
        mutations = variant.get("mutations")
        if not isinstance(mutations, list) or not mutations:
            raise FoldXError(f"Variant {variant_id}: mutations must be a nonempty list")
        sites: set[tuple] = set()
        normalized = []
        for mutation in mutations:
            _object(mutation, f"{variant_id} mutation", {"chain", "residue_number", "from", "to"})
            chain = mutation.get("chain")
            if not isinstance(chain, str) or chain not in mapping:
                raise FoldXError(f"Variant {variant_id}: mutation chain must be explicitly selected in chain_map")
            number = _integer(mutation.get("residue_number"), "residue_number (original PDB numbering; no renumbering)", 1, 9999)
            for field in ("from", "to"):
                aa = mutation.get(field)
                if not isinstance(aa, str) or aa not in CANONICAL_AA:
                    raise FoldXError(f"Variant {variant_id}: {field} must be one uppercase canonical amino acid")
            if mutation["from"] == mutation["to"]:
                raise FoldXError(f"Variant {variant_id}: identity substitutions are not mutation ΔΔG designs")
            if (chain, number) in sites:
                raise FoldXError(f"Variant {variant_id}: duplicate mutation site {chain}:{number}")
            sites.add((chain, number))
            normalized.append((chain, number, mutation["from"], mutation["to"]))
        design = tuple(sorted(normalized))
        if design in designs:
            raise FoldXError(f"Variant {variant_id}: duplicate mutation design")
        designs.add(design)
    foldx = _object(config.get("foldx", {}), "foldx", {"runs", "temperature", "pH", "ion_strength", "version", "timeout_seconds", "resources_dir"})
    version = foldx.get("version")
    if not isinstance(version, str) or not version.strip():
        raise FoldXError("foldx.version must declare the user-supplied FoldX version (not automatically verified)")
    normalized_foldx = {
        "runs": _integer(foldx.get("runs", 5), "foldx.runs", 1, 1000),
        "temperature": _number(foldx.get("temperature", 298), "foldx.temperature (K)", 1, 1000),
        "pH": _number(foldx.get("pH", 7), "foldx.pH", 0, 14),
        "ion_strength": _number(foldx.get("ion_strength", .05), "foldx.ion_strength (M)", 0, 5),
        "version": version.strip(),
        "timeout_seconds": _integer(foldx.get("timeout_seconds", 3600), "foldx.timeout_seconds", 1, 86400),
    }
    if "resources_dir" in foldx:
        value = foldx["resources_dir"]
        if not isinstance(value, (str, Path)) or not str(value).strip():
            raise FoldXError("foldx.resources_dir must be a directory path")
        resources = Path(value).expanduser().resolve()
        if not resources.is_dir():
            raise FoldXError(f"foldx.resources_dir not found: {resources}")
        normalized_foldx["resources_dir"] = str(resources)
    result["foldx"] = normalized_foldx
    return result


def _resources(config: dict) -> list[dict]:
    directory = config["foldx"].get("resources_dir")
    if not directory:
        return []
    root = Path(directory)
    inventory = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink() or (not path.is_file() and not path.is_dir()):
            raise FoldXError(f"Resources must contain regular files/directories, not links or special files: {path}")
        if path.is_file():
            relative = path.relative_to(root).as_posix()
            if relative.split("/")[0] in {"reference.pdb", "reference_Repair.pdb", "individual_list.txt", "stdout.log", "stderr.log"} or path.suffix.lower() in {".fxout", ".pdb"}:
                raise FoldXError(f"Resource name could collide with study input/output: {relative}")
            inventory.append({"path": relative, "sha256": _hash(path), "bytes": path.stat().st_size})
    return inventory


def _flags(foldx: dict) -> list[str]:
    return [f"--temperature={foldx['temperature']:g}", f"--pH={foldx['pH']:g}", f"--ionStrength={foldx['ion_strength']:g}"]


def plan_study(config: dict[str, Any], outdir: str | Path) -> dict[str, Any]:
    """Validate a point-mutation study and snapshot its inputs in a NEW folder."""
    config = validate_config(config)
    outdir = Path(outdir).expanduser().resolve()
    if outdir.exists():
        raise FoldXError(f"Output directory already exists; use a new directory to avoid stale results: {outdir}")
    resources = _resources(config)
    reference = config["reference"]
    fasta_source, pdb_source = Path(reference["fasta"]), Path(reference["pdb"])
    source_hashes = {"fasta_sha256": _hash(fasta_source), "pdb_sha256": _hash(pdb_source)}
    outdir.mkdir(parents=True, exist_ok=False)
    prepared = outdir / "prepared"
    prepared.mkdir()
    shutil.copyfile(fasta_source, prepared / "reference.fasta")
    shutil.copyfile(pdb_source, prepared / "source.pdb")
    if _hash(prepared / "reference.fasta") != source_hashes["fasta_sha256"] or _hash(prepared / "source.pdb") != source_hashes["pdb_sha256"]:
        raise FoldXError("Input changed while being copied; choose new output directory and retry")
    fasta = read_fasta(prepared / "reference.fasta")
    audit = prepare_pdb(prepared / "source.pdb", prepared / "reference.pdb", reference["chain_map"], fasta)
    variants = []
    for index, variant in enumerate(config["variants"], 1):
        tokens = []
        for mutation in sorted(variant["mutations"], key=lambda item: (item["chain"], item["residue_number"])):
            key = f"{mutation['chain']}:{mutation['residue_number']}"
            observed = audit["residue_map"].get(key)
            if observed != mutation["from"]:
                raise FoldXError(f"Variant {variant['id']}: WT mismatch at PDB {key}: expected {mutation['from']}, observed {observed!r}; check chain and original PDB numbering")
            tokens.append(f"{mutation['from']}{mutation['chain']}{mutation['residue_number']}{mutation['to']}")
        variants.append({"id": variant["id"], "foldx_index": index, "mutation_line": ",".join(tokens) + ";"})
    mutant_file = outdir / "individual_list.txt"
    mutant_file.write_text("\n".join(item["mutation_line"] for item in variants) + "\n", encoding="ascii")
    warnings = [
        "No structure prediction, missing-sidechain reconstruction or model-confidence assessment is performed; review the input structure before interpreting energies.",
        "Repeats sample FoldX algorithm behavior on one repaired backbone; seeds are not controlled and repeats are not guaranteed statistically independent.",
        "Positive mutant-minus-WT ΔΔG indicates destabilization within this model; no absolute stability or Tm ranking across different proteins/topologies is supported.",
    ]
    if config["synthetic"]:
        warnings.insert(0, "SYNTHETIC study: toy coordinates/fake outputs are software tests and have no scientific stability meaning.")
    plan = {
        "schema_version": 1, "study_id": config["study_id"], "comparison_type": config["comparison_type"], "synthetic": config["synthetic"],
        "config": config, "reference": {**source_hashes, "prepared_pdb_sha256": _hash(prepared / "reference.pdb"), "chain_map": reference["chain_map"]},
        "variants": variants, "mutation_file_sha256": _hash(mutant_file), "structure_audit": audit,
        "resources": resources, "warnings": warnings, "documentation": SOURCES,
        "commands": {
            "repair": ["<FOLDX_EXECUTABLE>", "--command=RepairPDB", "--pdb=reference.pdb", *_flags(config["foldx"])],
            "build_per_repeat": ["<FOLDX_EXECUTABLE>", "--command=BuildModel", "--pdb=reference_Repair.pdb", "--mutant-file=individual_list.txt", "--numberOfRuns=1", "--out-pdb=true", *_flags(config["foldx"])],
        },
    }
    _write_json(outdir / "plan.json", plan)
    return plan


def _copy_resources(config: dict, inventory: list[dict], destination: Path) -> None:
    if not inventory:
        return
    root = Path(config["foldx"]["resources_dir"])
    for item in inventory:
        source, target = root / item["path"], destination / item["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            raise FoldXError(f"Resource collision: {target}")
        shutil.copyfile(source, target)
        if _hash(target) != item["sha256"]:
            raise FoldXError(f"Resource changed since planning: {source}")


def _execute(command: list[str], cwd: Path, timeout: int, step: dict) -> None:
    step.update({"command": command, "cwd": str(cwd), "started_at": datetime.now(timezone.utc).isoformat()})
    stdout, stderr = cwd / "stdout.log", cwd / "stderr.log"
    try:
        with stdout.open("wb") as out, stderr.open("wb") as err:
            process = subprocess.run(command, cwd=cwd, stdout=out, stderr=err, timeout=timeout, check=False, shell=False, stdin=subprocess.DEVNULL)
        step["returncode"] = process.returncode
        if process.returncode != 0:
            raise FoldXError(f"FoldX exited with code {process.returncode} in {cwd}; see stdout.log and stderr.log")
    except subprocess.TimeoutExpired as exc:
        step["timed_out"] = True
        raise FoldXError(f"FoldX exceeded {timeout} seconds in {cwd}; see stdout.log and stderr.log") from exc
    finally:
        step["finished_at"] = datetime.now(timezone.utc).isoformat()
        step["files"] = [{"path": path.relative_to(cwd).as_posix(), "sha256": _hash(path), "bytes": path.stat().st_size} for path in sorted(cwd.rglob("*")) if path.is_file()]


def run_study(config: dict[str, Any], outdir: str | Path, executable: str | Path, license_confirmed: bool = False) -> dict[str, Any]:
    """Run one repair then isolated one-run BuildModel repeats.

    A user-confirmed valid license and supplied executable are necessary. No
    executable is downloaded and no original input/resource file is modified.
    The return manifest is written also on errors; errors are then re-raised.
    """
    plan = plan_study(config, outdir)
    root = Path(outdir).expanduser().resolve()
    config = plan["config"]
    manifest = {
        "schema_version": 1, "study_id": plan["study_id"], "comparison_type": plan["comparison_type"], "synthetic": plan["synthetic"],
        "status": "planned", "created_at": datetime.now(timezone.utc).isoformat(),
        "reference": copy.deepcopy(plan["reference"]), "variants": plan["variants"], "foldx": {**config["foldx"], "version_source": "user_declared", "seed_policy": "FoldX default; explicit seed not controlled"},
        "resources": plan["resources"], "warnings": plan["warnings"], "runs": [], "steps": [],
        "runtime": {"python": sys.version, "platform": platform.platform(), "machine": platform.machine(), "package": "protein_stability", "package_version": __version__, "source_sha256": {path.name: _hash(path) for path in sorted(Path(__file__).parent.glob("*.py"))}, "plan_sha256": _hash(root / "plan.json")},
    }
    manifest_path = root / "run_manifest.json"
    try:
        if license_confirmed is not True:
            raise FoldXError("Execution requires license_confirmed=True after obtaining FoldX and confirming the applicable license; the offline plan has been saved")
        if not isinstance(executable, (str, Path)) or not str(executable).strip():
            raise FoldXError("Supply the licensed FoldX executable path")
        executable_text = str(executable)
        resolved = shutil.which(executable_text) if "/" not in executable_text else executable_text
        binary = _file(resolved or executable_text, "FoldX executable")
        if not os.access(binary, os.X_OK):
            raise FoldXError(f"FoldX executable is not executable: {binary}")
        manifest["executable"] = {"path": str(binary), "sha256": _hash(binary), "bytes": binary.stat().st_size}
        manifest["license_confirmed"] = True
        manifest["status"] = "running"
        _write_json(manifest_path, manifest)
        repair_dir = root / "repair"
        repair_dir.mkdir()
        _copy_resources(config, plan["resources"], repair_dir)
        shutil.copyfile(root / "prepared/reference.pdb", repair_dir / "reference.pdb")
        step = {"name": "RepairPDB"}
        manifest["steps"].append(step)
        _execute([str(binary), *plan["commands"]["repair"][1:]], repair_dir, config["foldx"]["timeout_seconds"], step)
        repaired = _file(repair_dir / "reference_Repair.pdb", "RepairPDB repaired structure")
        # Revalidate full sequence/numbering after the external tool; never accept
        # a repaired file silently truncated or changed at the mutation sites.
        repaired_audit = prepare_pdb(repaired, root / "prepared/repaired_validated.pdb", config["reference"]["chain_map"], read_fasta(root / "prepared/reference.fasta"))
        if repaired_audit["residue_map"] != plan["structure_audit"]["residue_map"]:
            raise FoldXError("RepairPDB changed residue identities/numbering; aborting")
        repaired_hash = _hash(repaired)
        manifest["reference"]["repaired_pdb_sha256"] = repaired_hash
        replicas = root / "replicates"
        replicas.mkdir()
        for replicate in range(1, config["foldx"]["runs"] + 1):
            work = replicas / f"run_{replicate:03d}"
            work.mkdir()
            _copy_resources(config, plan["resources"], work)
            shutil.copyfile(repaired, work / "reference_Repair.pdb")
            if _hash(work / "reference_Repair.pdb") != repaired_hash:
                raise FoldXError("Repaired reference changed before BuildModel; aborting")
            shutil.copyfile(root / "individual_list.txt", work / "individual_list.txt")
            if _hash(work / "individual_list.txt") != plan["mutation_file_sha256"]:
                raise FoldXError("Mutation file changed since planning; aborting")
            run = {"replicate": replicate, "directory": work.relative_to(root).as_posix(), "variants": [{"variant_id": item["id"], "foldx_index": item["foldx_index"], "expected_model_id": f"reference_Repair_{item['foldx_index']}.pdb"} for item in plan["variants"]]}
            run["mutation_file_sha256"] = plan["mutation_file_sha256"]
            run["repaired_pdb_sha256"] = repaired_hash
            manifest["runs"].append(run)
            step = {"name": "BuildModel", "replicate": replicate}
            manifest["steps"].append(step)
            _execute([str(binary), *plan["commands"]["build_per_repeat"][1:]], work, config["foldx"]["timeout_seconds"], step)
            differences = sorted(path for path in work.glob("*.fxout") if path.name.startswith(("Dif_", "Diff_")))
            if len(differences) != 1:
                raise FoldXError(f"Expected exactly one Dif_/Diff_ output in {work}; found {len(differences)}")
            for prefix in ("Raw_", "Average_", "PdbList_"):
                matches = list(work.glob(prefix + "*.fxout"))
                if len(matches) != 1 or matches[0].stat().st_size == 0:
                    raise FoldXError(f"Incomplete BuildModel output in {work}: expected one nonempty {prefix}*.fxout")
            rows = parse_diff(differences[0])
            expected = [item["expected_model_id"].removesuffix(".pdb") for item in run["variants"]]
            observed = [row["model_id"].removesuffix(".pdb") for row in rows]
            if observed != expected:
                raise FoldXError(f"BuildModel output row/model mapping mismatch in {work}: expected {expected}, observed {observed}")
            for item in run["variants"]:
                model = work / item["expected_model_id"]
                wt = work / ("WT_" + item["expected_model_id"])
                _file(model, "BuildModel mutant structure")
                _file(wt, "BuildModel paired WT structure")
                expected_mutant = dict(plan["structure_audit"]["residue_map"])
                for mutation in config["variants"][item["foldx_index"] - 1]["mutations"]:
                    expected_mutant[f"{mutation['chain']}:{mutation['residue_number']}"] = mutation["to"]
                if inspect_pdb(model)["residue_map"] != expected_mutant:
                    raise FoldXError(f"BuildModel mutant identity/numbering mismatch: {model}")
                if inspect_pdb(wt)["residue_map"] != plan["structure_audit"]["residue_map"]:
                    raise FoldXError(f"BuildModel paired WT identity/numbering mismatch: {wt}")
            run["diff_file"] = differences[0].relative_to(root).as_posix()
            run["diff_sha256"] = _hash(differences[0])
            _write_json(manifest_path, manifest)
        manifest["status"] = "completed"
        manifest["finished_at"] = datetime.now(timezone.utc).isoformat()
        _write_json(manifest_path, manifest)
        return manifest
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        manifest["finished_at"] = datetime.now(timezone.utc).isoformat()
        _write_json(manifest_path, manifest)
        raise
