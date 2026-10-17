import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from protein_prediction.pipeline import fasta, pdb, prepare, parse_dif, summarize

ROOT = Path(__file__).resolve().parents[1]

class WorkflowTests(unittest.TestCase):
    def test_demo_and_report(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'job'
            m=prepare(ROOT/'examples/config.json',out)
            self.assertEqual(m['status'],'planned')
            self.assertIn('--numberOfRuns=3',m['commands'][1])
            rows=summarize(ROOT/'examples/Dif_synthetic.fxout',out/'manifest.json',out/'report')
            self.assertAlmostEqual(rows[0]['mean_ddg'],1.2)
            self.assertAlmostEqual(rows[0]['sd'],.2)
            self.assertIn('synthetic', (out/'report/index.html').read_text())
            with self.assertRaises(FileExistsError):
                prepare(ROOT/'examples/config.json',out)

    def test_fasta_rejections(self):
        for content in ['AA','>a\nAX','>a\nAA\n>a\nAC','>a\n']:
            with self.subTest(content=content), tempfile.TemporaryDirectory() as d:
                p=Path(d)/'x';p.write_text(content)
                with self.assertRaises(ValueError):fasta(p)

    def test_mutation_and_comparison_gates(self):
        base=json.loads((ROOT/'examples/config.json').read_text())
        base['pdb']=str(ROOT/'examples/toy.pdb');base['fasta']=str(ROOT/'examples/toy.fasta')
        for change in [{'comparison':'linker'},{'repeats':1},{'chains':{'A':'wrong'}}, {'variants':[{'id':'bad','mutations':['CA1V']}]}, {'variants':[{'id':'bad','mutations':['AA1V','AA1L']}]}]:
            with self.subTest(change=change), tempfile.TemporaryDirectory() as d:
                config=Path(d)/'c.json';config.write_text(json.dumps(base|change))
                with self.assertRaises(ValueError):prepare(config,Path(d)/'job')

    def test_pdb_gate(self):
        text=(ROOT/'examples/toy.pdb').read_text()
        for value in [text+'HETATM\n', text.replace('ATOM  ', 'HETATM',1), text.replace(' O   ALA',' XX  ALA'), text+text]:
            with self.subTest(value=value[:20]),tempfile.TemporaryDirectory() as d:
                p=Path(d)/'bad.pdb';p.write_text(value)
                with self.assertRaises(ValueError):pdb(p)

    def test_fxout_gate(self):
        for name,content in [('Raw_x.fxout',''),('Dif_x.fxout','Pdb\tTotal Energy\na\tnan\n'),('Dif_x.fxout','Pdb\tTotal Energy\na\t1\t2\n')]:
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/name;p.write_text(content)
                with self.assertRaises(ValueError):parse_dif(p)

    def test_incomplete_duplicate_and_unknown_rows(self):
        for rows in ['x_1_1.pdb\t1\n','x_1_1.pdb\t1\nx_1_1.pdb\t2\n','unknown\t1\n']:
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/'Dif_x.fxout';p.write_text('Pdb\tTotal Energy\n'+rows)
                manifest=Path(d)/'manifest.json';manifest.write_text(json.dumps({'config':json.loads((ROOT/'examples/config.json').read_text())}))
                with self.assertRaises(ValueError):summarize(p,manifest,Path(d)/'report')

    def test_missing_binary_records_failure(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'job'
            with self.assertRaises(ValueError):prepare(ROOT/'examples/config.json',out,True,'/nonexistent/foldx')
            self.assertEqual(json.loads((out/'manifest.json').read_text())['status'],'failed')

    def test_mock_execution_contract(self):
        # Integration contract only: simulated executable output is never scientific evidence.
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'job'
            def fake_run(command,**kwargs):
                if '--command=RepairPDB' in command:
                    (out/'input_Repair.pdb').write_bytes((out/'input.pdb').read_bytes())
                if '--command=BuildModel' in command:
                    (out/'Dif_input_Repair.fxout').write_bytes((ROOT/'examples/Dif_synthetic.fxout').read_bytes())
                import subprocess
                return subprocess.CompletedProcess(command,0,stdout='mock help',stderr='')
            with patch('protein_prediction.pipeline.shutil.which',return_value='/usr/bin/true'),patch('protein_prediction.pipeline.subprocess.run',side_effect=fake_run):
                m=prepare(ROOT/'examples/config.json',out,True)
            self.assertEqual(m['status'],'completed')
            self.assertTrue((out/'report/index.html').exists())

if __name__=='__main__':unittest.main()
