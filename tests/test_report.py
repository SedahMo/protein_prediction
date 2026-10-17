import csv
import json
from pathlib import Path
import tempfile
import unittest

from protein_stability.report import blocked_report, make_demo, write_report
from protein_stability.results import collect_results


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_demo_collects_two_variants_five_repeats(self):
        paths = make_demo(self.root)
        self.assertEqual(set(paths), {"results_json", "summary_csv", "report_html", "run_manifest"})
        data = json.loads(Path(paths["results_json"]).read_text())
        self.assertTrue(data["synthetic"])
        self.assertEqual(len(data["variants"]), 2)
        self.assertTrue(all(v["n"] == 5 for v in data["variants"]))
        self.assertAlmostEqual(data["variants"][0]["mean_ddg_kcal_mol"], -0.3)
        self.assertAlmostEqual(data["variants"][1]["mean_ddg_kcal_mol"], 0.4)
        self.assertEqual(collect_results(self.root / "synthetic_run"), data)
        document = Path(paths["report_html"]).read_text()
        self.assertIn("SYNTHETIC DEMONSTRATION", document)
        self.assertIn("No biological stability calculation", document)
        self.assertIn("<svg", document)
        self.assertNotIn("<script", document)
        self.assertNotIn('src="https://', document)
        self.assertIn("not experimental uncertainty", document)
        with Path(paths["summary_csv"]).open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["synthetic"] == "True" for row in rows))

    def test_blocked_and_nonpoint_templates_do_not_plot(self):
        for mode in ("point_mutation", "linker", "domain_order", "cross_protein"):
            with self.subTest(mode=mode):
                data = blocked_report("study", "Supply validated structural models.", mode)
                paths = write_report(data, self.root / mode)
                document = Path(paths["report_html"]).read_text()
                self.assertNotIn("<svg", document)
                self.assertIn("ANALYSIS TEMPLATE", document)
                self.assertIn("No completed energy calculations", document)
                self.assertIn("Supply validated structural models.", document)
                self.assertIn("Comparison rules", document)
                self.assertIn("cross_protein", document)
                self.assertEqual(json.loads(Path(paths["results_json"]).read_text())["variants"], [])

    def test_nonpoint_energy_summary_cannot_get_a_ranking(self):
        paths = make_demo(self.root / "source")
        data = json.loads(Path(paths["results_json"]).read_text())
        data["comparison_type"] = "linker"
        with self.assertRaisesRegex(ValueError, "not eligible"):
            write_report(data, self.root / "bad")
        self.assertFalse((self.root / "bad").exists())

    def test_html_escapes_input_and_csv_avoids_formulas(self):
        paths = make_demo(self.root / "source")
        data = json.loads(Path(paths["results_json"]).read_text())
        data["study_id"] = '<script>alert("x")</script>'
        data["variants"][0]["variant_id"] = "=1+1"
        paths = write_report(data, self.root / "escaped")
        self.assertNotIn('<script>alert', Path(paths["report_html"]).read_text())
        self.assertIn("&lt;script&gt;", Path(paths["report_html"]).read_text())
        with Path(paths["summary_csv"]).open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]["variant_id"], "'=1+1")

    def test_nonfinite_json_data_is_rejected_before_write(self):
        data = blocked_report("study", "reason")
        data["invalid"] = float("nan")
        with self.assertRaises(ValueError):
            write_report(data, self.root / "bad")
        self.assertFalse((self.root / "bad").exists())


if __name__ == "__main__":
    unittest.main()
