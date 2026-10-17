import hashlib
import json
from pathlib import Path
import statistics
import tempfile
import unittest

from protein_stability.results import ResultsError, collect_results, parse_diff


class ParseDiffTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def parse(self, content, filename="Dif_reference_Repair.fxout"):
        path = self.root / filename
        path.write_bytes(content.encode("utf-8"))
        return parse_diff(path)

    def test_realistic_preamble_bom_crlf_terms(self):
        rows = self.parse("\ufeffFoldX 5 output\r\nBuildModel\r\n\r\nPdb\ttotal energy\tBackbone Hbond\tSidechain Hbond\r\nreference_Repair_1\t-1.25\t0.2\t-0.7\r\n\r\n")
        self.assertEqual(rows, [{"model_id": "reference_Repair_1", "ddg_kcal_mol": -1.25, "terms": {"Backbone Hbond": 0.2, "Sidechain Hbond": -0.7}}])

    def test_both_official_and_alternate_difference_prefix(self):
        self.assertEqual(self.parse("Pdb\ttotal energy\nm1\t0.5\n", "Diff_reference.fxout")[0]["ddg_kcal_mol"], 0.5)

    def test_consistent_single_terminal_tab(self):
        self.assertEqual(self.parse("Pdb\ttotal energy\tterm\t\nm1\t1\t2\t\n")[0]["terms"], {"term": 2.0})

    def test_never_parse_raw_or_average_as_ddg(self):
        for name in ("Raw_reference.fxout", "Average_reference.fxout", "energy.txt"):
            with self.subTest(name=name), self.assertRaises(ResultsError):
                self.parse("Pdb\ttotal energy\nm1\t0.5\n", name)

    def test_missing_or_wrong_header(self):
        for content in ("", "Pdb total energy\nm1 0.5\n", "Pdb\tother energy\nm1\t0.5\n", "Pdb\ttotal energy\n"):
            with self.subTest(content=content), self.assertRaises(ResultsError):
                self.parse(content)

    def test_duplicate_headers_repeated_header_duplicate_model(self):
        for content in ("Pdb\ttotal energy\tTotal Energy\nm1\t1\t2\n", "Pdb\ttotal energy\nm1\t1\nm1\t2\n", "Pdb\ttotal energy\nPdb\ttotal energy\n"):
            with self.subTest(content=content), self.assertRaises(ResultsError):
                self.parse(content)

    def test_truncated_extra_empty_nonnumeric_nonfinite(self):
        for row in ("m1\t1", "m1\t1\t2\textra", "m1\t1\t", "m1\tbad\t2", "m1\tnan\t2", "m1\t1\t-inf", "m1\t1e999\t2", "m1\t1_000\t2"):
            with self.subTest(row=row), self.assertRaises(ResultsError):
                self.parse("Pdb\ttotal energy\tterm\n" + row + "\n")

    def test_inconsistent_or_multiple_terminal_tabs_rejected(self):
        for content in ("Pdb\ttotal energy\t\nm1\t1\n", "Pdb\ttotal energy\t\t\nm1\t1\t\t\n", "Pdb\ttotal energy\nm1\t1\t\n"):
            with self.subTest(content=content), self.assertRaises(ResultsError):
                self.parse(content)


class CollectResultsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        digest = "a" * 64
        self.manifest = {"schema_version": 1, "study_id": "study-1", "comparison_type": "point_mutation", "synthetic": True, "status": "completed", "reference": {"pdb_sha256": digest, "fasta_sha256": digest, "prepared_pdb_sha256": digest}, "foldx": {"runs": 2}, "variants": [{"id": "A1V", "foldx_index": 1, "mutation_line": "AA1V;"}, {"id": "G2S", "foldx_index": 2, "mutation_line": "GA2S;"}], "runs": []}
        for replicate, values in ((1, (-1.0, 2.0)), (2, (-2.0, 4.0))):
            directory = f"replicates/run_{replicate:03d}"
            (self.root / directory).mkdir(parents=True)
            diff_file = f"{directory}/Dif_reference_Repair.fxout"
            path = self.root / diff_file
            path.write_text("FoldX fixture\nPdb\ttotal energy\tterm\n" + "".join(f"reference_Repair_{index}\t{value}\t0\n" for index, value in enumerate(values, 1)), encoding="utf-8")
            self.manifest["runs"].append({"replicate": replicate, "directory": directory, "diff_file": diff_file, "diff_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "variants": [{"variant_id": v["id"], "foldx_index": v["foldx_index"], "expected_model_id": f"reference_Repair_{v['foldx_index']}.pdb"} for v in self.manifest["variants"]]})
        self.save()

    def save(self):
        (self.root / "run_manifest.json").write_text(json.dumps(self.manifest), encoding="utf-8")

    def collect(self):
        self.save()
        return collect_results(self.root)

    def genuine(self):
        self.manifest["synthetic"] = False
        for relative, key in (("prepared/source.pdb", "pdb_sha256"), ("prepared/reference.fasta", "fasta_sha256"), ("prepared/reference.pdb", "prepared_pdb_sha256"), ("repair/reference_Repair.pdb", "repaired_pdb_sha256")):
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("test snapshot for " + relative, encoding="utf-8")
            self.manifest["reference"][key] = hashlib.sha256(path.read_bytes()).hexdigest()
        mutation_path = self.root / "individual_list.txt"
        mutation_path.write_text("AA1V;\nGA2S;\n", encoding="utf-8")
        plan = {"schema_version": 1, "study_id": self.manifest["study_id"], "comparison_type": self.manifest["comparison_type"], "synthetic": False, "variants": self.manifest["variants"], "reference": {k: v for k, v in self.manifest["reference"].items() if k != "repaired_pdb_sha256"}, "mutation_file_sha256": hashlib.sha256(mutation_path.read_bytes()).hexdigest(), "config": {"foldx": dict(self.manifest["foldx"])}}
        plan_path = self.root / "plan.json"
        plan_path.write_text(json.dumps(plan), encoding="utf-8")
        self.manifest["runtime"] = {"plan_sha256": hashlib.sha256(plan_path.read_bytes()).hexdigest()}
        return plan

    def test_genuine_run_validates_local_snapshots_and_plan(self):
        self.genuine()
        self.assertFalse(self.collect()["synthetic"])

    def test_genuine_run_requires_complete_folder(self):
        self.genuine()
        (self.root / "repair/reference_Repair.pdb").unlink()
        with self.assertRaisesRegex(ResultsError, "Complete run folder"):
            self.collect()

    def test_genuine_modified_prepared_snapshot_rejected(self):
        self.genuine()
        (self.root / "prepared/reference.pdb").write_text("changed", encoding="utf-8")
        with self.assertRaisesRegex(ResultsError, "SHA-256 mismatch"):
            self.collect()

    def test_genuine_plan_hash_and_scope_rejected(self):
        plan = self.genuine()
        plan["study_id"] = "different-study"
        path = self.root / "plan.json"
        path.write_text(json.dumps(plan), encoding="utf-8")
        with self.assertRaisesRegex(ResultsError, "SHA-256 mismatch"):
            self.collect()
        self.manifest["runtime"]["plan_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        with self.assertRaisesRegex(ResultsError, "disagree on study_id"):
            self.collect()

    def test_recorded_step_output_hash_rejected_after_modification(self):
        run = self.manifest["runs"][0]
        path = self.root / run["directory"] / "model.pdb"
        path.write_text("before", encoding="utf-8")
        self.manifest["steps"] = [{"name": "BuildModel", "replicate": 1, "returncode": 0, "cwd": "/old/archive/location", "files": [{"path": "model.pdb", "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size}]}]
        self.collect()
        path.write_text("after", encoding="utf-8")
        with self.assertRaisesRegex(ResultsError, "SHA-256 mismatch"):
            self.collect()

    def test_complete_scope_and_sample_sd(self):
        data = self.collect()
        first, second = data["variants"]
        self.assertEqual(first["variant_id"], "A1V")
        self.assertEqual(first["mean_ddg_kcal_mol"], -1.5)
        self.assertAlmostEqual(first["sd_ddg_kcal_mol"], statistics.stdev([-1, -2]))
        self.assertEqual((first["min_ddg_kcal_mol"], first["max_ddg_kcal_mol"], first["n"]), (-2, -1, 2))
        self.assertEqual(second["mean_ddg_kcal_mol"], 3.0)
        self.assertEqual(len(data["provenance"]["source_files"]), 2)
        self.assertIn("not experimental", data["uncertainty_label"])
        self.assertIn("positive", data["sign_convention"])

    def test_single_repeat_sd_is_null(self):
        self.manifest["foldx"]["runs"] = 1
        self.manifest["runs"] = self.manifest["runs"][:1]
        self.assertIsNone(self.collect()["variants"][0]["sd_ddg_kcal_mol"])

    def test_incomplete_replicates_and_status(self):
        self.manifest["runs"].pop()
        with self.assertRaises(ResultsError):
            self.collect()
        self.manifest["status"] = "running"
        with self.assertRaises(ResultsError):
            self.collect()

    def test_output_hash_change_is_rejected(self):
        path = self.root / self.manifest["runs"][0]["diff_file"]
        path.write_text(path.read_text().replace("-1.0", "-3.0"), encoding="utf-8")
        with self.assertRaisesRegex(ResultsError, "SHA-256"):
            self.collect()

    def test_missing_or_extra_rows_rejected(self):
        run = self.manifest["runs"][0]
        path = self.root / run["diff_file"]
        path.write_text("Pdb\ttotal energy\nreference_Repair_1\t1\n", encoding="utf-8")
        run["diff_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        with self.assertRaisesRegex(ResultsError, "expected 2"):
            self.collect()

    def test_row_order_is_not_silently_reassigned(self):
        run = self.manifest["runs"][0]
        path = self.root / run["diff_file"]
        path.write_text("Pdb\ttotal energy\nreference_Repair_2\t2\nreference_Repair_1\t1\n", encoding="utf-8")
        run["diff_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        with self.assertRaisesRegex(ResultsError, "does not match"):
            self.collect()

    def test_inconsistent_mapping_rejected(self):
        self.manifest["runs"][0]["variants"][0]["variant_id"] = "wrong_variant"
        with self.assertRaises(ResultsError):
            self.collect()

    def test_output_hash_must_be_recorded(self):
        self.manifest["runs"][0].pop("diff_sha256")
        with self.assertRaisesRegex(ResultsError, "diff_sha256"):
            self.collect()

    def test_cross_backbone_or_study_pooling_rejected(self):
        for key, value in (("reference", {"pdb_sha256": "b" * 64}), ("study_id", "another-study"), ("comparison_type", "linker")):
            with self.subTest(key=key):
                run = self.manifest["runs"][0]
                run[key] = value
                with self.assertRaisesRegex(ResultsError, "scope differs"):
                    self.collect()
                del run[key]

    def test_nonpoint_comparison_rejected(self):
        for mode in ("linker", "domain_order", "cross_protein"):
            self.manifest["comparison_type"] = mode
            with self.subTest(mode=mode), self.assertRaises(ResultsError):
                self.collect()

    def test_duplicate_replica_or_output_rejected(self):
        self.manifest["runs"][1]["replicate"] = 1
        with self.assertRaises(ResultsError):
            self.collect()
        self.manifest["runs"][1]["replicate"] = 2
        self.manifest["runs"][1]["directory"] = self.manifest["runs"][0]["directory"]
        self.manifest["runs"][1]["diff_file"] = self.manifest["runs"][0]["diff_file"]
        with self.assertRaises(ResultsError):
            self.collect()

    def test_path_traversal_rejected(self):
        self.manifest["runs"][0]["diff_file"] = "../Dif_outside.fxout"
        with self.assertRaisesRegex(ResultsError, "escapes"):
            self.collect()

    def test_duplicate_manifest_json_key_rejected(self):
        (self.root / "run_manifest.json").write_text('{"schema_version":1,"schema_version":1}', encoding="utf-8")
        with self.assertRaisesRegex(ResultsError, "Duplicate JSON"):
            collect_results(self.root)


if __name__ == "__main__":
    unittest.main()
