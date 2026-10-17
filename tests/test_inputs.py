"""Input tests use synthetic PDB coordinates with meaningless geometry.

Fixtures exercise text format and identity checks only; they must never be used
as structure-prediction or energetic-validation data.
"""

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from protein_stability.inputs import inspect_fasta, inspect_pdb, prepare_pdb, read_fasta


def atom(serial, name, residue="GLY", chain="A", number=1, *, altloc=" ",
         insertion=" ", x=1.0, record="ATOM"):
    """Format one synthetic atom; coordinates have no structural meaning."""
    element = name[0]
    return (
        f"{record:<6}{serial:5d} {name:>4}{altloc}{residue:>3} {chain}{number:4d}{insertion}   "
        f"{x:8.3f}{2.0:8.3f}{3.0:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2}\n"
    )


def backbone(residue="GLY", chain="A", number=1, *, start=1, insertion=" "):
    return "".join(atom(start + i, name, residue, chain, number, insertion=insertion)
                   for i, name in enumerate(("N", "CA", "C", "O")))


class InputTestCase(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, text, name="input.pdb"):
        path = self.root / name
        path.write_text(text, encoding="ascii")
        return path


class FastaTests(InputTestCase):
    def test_wrapped_lowercase_and_descriptive_headers(self):
        path = self.write(">first description\nacd\nefghik\n\n>second\nLMNPQRSTVWY\n", "input.fa")
        self.assertEqual(read_fasta(path), {"first": "ACDEFGHIK", "second": "LMNPQRSTVWY"})

    def test_rejects_invalid_records(self):
        invalid = ("", "ACD\n", ">\nACD\n", ">one\n", ">one\n>two\nACD\n",
                   ">one\nACD\n>one duplicate\nACD\n", ">one\nAC-D\n",
                   ">one\nACD*\n", ">one\nACDX\n", ">one\nACD1\n",
                   ">one\nAC D\n", ">one\nACDB\n", ">one\nACDU\n")
        for text in invalid:
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_fasta(self.write(text, "input.fa"))

    def test_inspection_is_json_serializable_without_stability_score(self):
        path = self.write(">scfv\ncggggsggggsggggs\n", "input.fa")
        report = inspect_fasta(path)
        json.dumps(report, allow_nan=False)
        self.assertEqual(report["record_count"], 1)
        record = report["records"][0]
        self.assertEqual(record["length"], 16)
        self.assertEqual(record["cysteine_count"], 1)
        self.assertEqual(record["sha256"], hashlib.sha256(b"CGGGGSGGGGSGGGGS").hexdigest())
        self.assertGreaterEqual(len(record["warnings"]), 2)
        self.assertNotIn("stability_score", record)


class PdbTests(InputTestCase):
    def test_selected_chain_identity_and_preserved_numbering(self):
        source = self.write(
            "HEADER    SYNTHETIC INPUT TEST; GEOMETRY IS MEANINGLESS\n" +
            backbone("GLY", "A", -1) + backbone("CYS", "A", 3, start=5) +
            backbone("ALA", "B", 7, start=9) +
            atom(13, "O", "HOH", "W", 1, record="HETATM") +
            atom(14, "C", "LIG", "L", 1, record="HETATM") + "END\n")
        destination = self.root / "prepared" / "model.pdb"
        audit = prepare_pdb(source, destination, {"A": "target"}, {"target": "GC"})
        self.assertEqual(audit["residue_map"], {"A:-1": "G", "A:3": "C"})
        self.assertEqual(audit["chains"]["A"]["residue_numbers"], [-1, 3])
        self.assertEqual(audit["chains"]["A"]["fasta_id"], "target")
        self.assertEqual(audit["removed"]["water_atoms"], 1)
        self.assertEqual(audit["removed"]["ligand_atoms"], 1)
        self.assertEqual(audit["removed"]["unselected_protein_atoms"], 4)
        self.assertTrue(any("discontinuity" in item for item in audit["warnings"]))
        output = destination.read_text()
        self.assertTrue(output.endswith("TER       9      CYS A   3 \nEND\n"))
        original_atom_lines = [line for line in source.read_text().splitlines()
                               if line.startswith("ATOM") and line[21] == "A"]
        self.assertEqual([line for line in output.splitlines() if line.startswith("ATOM")],
                         [line.rstrip() for line in original_atom_lines])
        self.assertEqual(audit["output_sha256"], hashlib.sha256(destination.read_bytes()).hexdigest())
        json.dumps(audit, allow_nan=False)

    def test_exact_sequence_match_required(self):
        source = self.write(backbone())
        for expected in ("A", "GG"):
            with self.subTest(expected=expected), self.assertRaisesRegex(ValueError, "exactly match"):
                prepare_pdb(source, self.root / "out.pdb", {"A": "target"}, {"target": expected})
        self.assertFalse((self.root / "out.pdb").exists())

    def test_requires_explicit_valid_chain_and_fasta_mapping(self):
        source = self.write(backbone())
        cases = (({}, {"target": "G"}), ({"B": "target"}, {"target": "G"}),
                 ({"A": "missing"}, {"target": "G"}), ({" ": "target"}, {"target": "G"}),
                 ({"AA": "target"}, {"target": "G"}), ({"A": "target"}, {"target": "g"}))
        for chain_map, fasta in cases:
            with self.subTest(chain_map=chain_map, fasta=fasta), self.assertRaises(ValueError):
                prepare_pdb(source, self.root / "out.pdb", chain_map, fasta)

    def test_refuses_source_overwrite(self):
        source = self.write(backbone())
        original = source.read_bytes()
        with self.assertRaisesRegex(ValueError, "differ from"):
            prepare_pdb(source, source, {"A": "target"}, {"target": "G"})
        self.assertEqual(source.read_bytes(), original)

    def test_single_explicit_model_accepted_and_ensemble_rejected(self):
        valid = self.write("MODEL        1\n" + backbone() + "ENDMDL\nEND\n")
        self.assertEqual(inspect_pdb(valid)["model_count"], 1)
        for text in ("MODEL        1\n" + backbone() + "ENDMDL\nMODEL        2\n" + backbone() + "ENDMDL\n",
                     "MODEL        1\n" + backbone(),
                     "ENDMDL\n" + backbone(),
                     backbone() + "MODEL        1\n" + backbone() + "ENDMDL\n",
                     "MODEL        1\n" + backbone() + "ENDMDL\n" + backbone()):
            with self.subTest(text=text), self.assertRaises(ValueError):
                inspect_pdb(self.write(text))

    def test_insertion_code_inspectable_but_not_executable(self):
        source = self.write(backbone(insertion="A"))
        self.assertEqual(inspect_pdb(source)["residue_map"], {"A:1A": "G"})
        with self.assertRaisesRegex(ValueError, "insertion codes"):
            prepare_pdb(source, self.root / "out.pdb", {"A": "target"}, {"target": "G"})

    def test_noncanonical_and_missing_backbone_rejected(self):
        for text in (backbone("MSE"), atom(1, "CA"), backbone().replace("GLY", "UNK")):
            with self.subTest(text=text), self.assertRaises(ValueError):
                inspect_pdb(self.write(text))

    def test_nonfinite_and_truncated_coordinates_rejected(self):
        for text in (atom(1, "N", x=float("nan")) + backbone(start=2),
                     atom(1, "N", x=float("inf")) + backbone(start=2),
                     "ATOM      1\n", "data_structure\n_atom_site.id\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                inspect_pdb(self.write(text))

    def test_duplicate_atoms_and_residue_identifiers_rejected(self):
        for text in (backbone() + atom(5, "CA"),
                     backbone(number=1) + backbone(number=2, start=5) + backbone(number=1, start=9),
                     backbone("GLY") + atom(5, "CB", "ALA")):
            with self.subTest(text=text), self.assertRaises(ValueError):
                inspect_pdb(self.write(text))

    def test_deterministic_altloc_selection_and_normalized_output(self):
        source = self.write(backbone("ALA") +
                            atom(5, "CB", "ALA", altloc="A", x=4.0) +
                            atom(6, "CB", "ALA", altloc="B", x=5.0))
        output = self.root / "out.pdb"
        audit = prepare_pdb(source, output, {"A": "target"}, {"target": "A"})
        self.assertEqual(audit["removed"]["alternate_location_atoms"], 1)
        selected = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]
        self.assertEqual(len(selected), 5)
        self.assertTrue(all(line[16] == " " for line in selected))
        self.assertEqual(float(selected[-1][30:38]), 4.0)

    def test_ambiguous_altloc_and_only_non_a_rejected(self):
        for tail in (atom(5, "CB", "ALA", altloc="B"),
                     atom(5, "CB", "ALA") + atom(6, "CB", "ALA", altloc="A"),
                     atom(5, "CB", "ALA", altloc="A") + atom(6, "CB", "ALA", altloc="A")):
            with self.subTest(tail=tail), self.assertRaises(ValueError):
                inspect_pdb(self.write(backbone("ALA") + tail))

    def test_unselected_noncanonical_chain_does_not_block_explicit_selection(self):
        source = self.write(backbone("GLY", "A") + backbone("MSE", "B", start=5))
        audit = prepare_pdb(source, self.root / "out.pdb", {"A": "target"}, {"target": "G"})
        self.assertEqual(audit["selected_chains"], ["A"])
        with self.assertRaisesRegex(ValueError, "noncanonical"):
            inspect_pdb(source)

    def test_ter_prevents_merging_disconnected_segments(self):
        for ter in ("TER\n", "TER       5      GLY A   1\n"):
            with self.subTest(ter=ter):
                source = self.write(backbone("GLY", "A", 1) + ter +
                                    backbone("ALA", "A", 2, start=6))
                with self.assertRaisesRegex(ValueError, "reappears after TER"):
                    prepare_pdb(source, self.root / "out.pdb", {"A": "target"}, {"target": "GA"})

    def test_ter_between_distinct_chains_remains_valid(self):
        source = self.write(backbone("GLY", "A", 1) + "TER\n" +
                            backbone("ALA", "B", 1, start=6) + "TER\n")
        report = prepare_pdb(source, self.root / "out.pdb", {"A": "first", "B": "second"},
                             {"first": "G", "second": "A"})
        self.assertEqual(report["selected_chains"], ["A", "B"])

    def test_short_ter_uses_latest_unselected_chain(self):
        source = self.write(backbone("GLY", "A", 1) + backbone("ALA", "B", 1, start=5) +
                            "TER\n" + backbone("CYS", "A", 2, start=10))
        report = prepare_pdb(source, self.root / "out.pdb", {"A": "target"}, {"target": "GC"})
        self.assertEqual(report["chains"]["A"]["sequence"], "GC")


if __name__ == "__main__":
    unittest.main()
