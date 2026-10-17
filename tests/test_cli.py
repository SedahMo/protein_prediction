import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from protein_stability.cli import doctor, load_config, main


class CliTests(unittest.TestCase):
    def test_config_paths_relative_to_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "study.json"
            path.write_text(json.dumps({"reference": {"fasta": "input.fa", "pdb": "model.pdb"},
                                        "foldx": {"resources_dir": "resources"}}))
            config = load_config(path)
            self.assertEqual(config["reference"]["fasta"], str(Path(tmp) / "input.fa"))
            self.assertEqual(config["foldx"]["resources_dir"], str(Path(tmp) / "resources"))

    def test_doctor_does_not_claim_execution_or_license(self):
        result = doctor("no_such_foldx_program_28913")
        self.assertTrue(result["offline_modules_ready"])
        self.assertFalse(result["foldx_detected"])
        self.assertFalse(result["foldx_execution_verified"])
        self.assertFalse(result["foldx_license_verified"])

    def test_duplicate_config_keys_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "study.json"
            path.write_text('{"study_id":"one","study_id":"two"}')
            with self.assertRaisesRegex(ValueError, "Duplicate JSON key"):
                load_config(path)

    def test_bad_fasta_cli_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.fasta"
            path.write_text(">bad\nAC*D\n")
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                self.assertEqual(main(["validate-fasta", str(path)]), 2)
            self.assertIn("ERROR:", err.getvalue())

    def test_comparison_template_has_no_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "report"
            with contextlib.redirect_stdout(io.StringIO()):
                status = main(["comparison-template", "--type", "domain_order", "--out", str(out)])
            self.assertEqual(status, 0)
            data = json.loads((out / "results.json").read_text())
            self.assertEqual(data["variants"], [])
            self.assertEqual(data["comparison_type"], "domain_order")


if __name__ == "__main__":
    unittest.main()
