"""Synthetic external-process fixtures test orchestration, never FoldX energies."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from protein_stability.foldx import plan_study, run_study, validate_config
from protein_stability.results import collect_results

FAKE = r'''
import pathlib, shutil, sys, time
mode = MODE
args = dict(arg[2:].split('=',1) for arg in sys.argv[1:])
print('SYNTHETIC FAKE PROCESS — NOT FOLDX')
print('synthetic stderr', file=sys.stderr)
if mode == 'exit': sys.exit(7)
if mode == 'timeout': time.sleep(3)
if args['command'] == 'RepairPDB':
    if mode != 'missing_repaired':
        data = pathlib.Path(args['pdb']).read_text()
        if mode == 'changed_repair': data = data.replace('ALA A   1', 'VAL A   1')
        pathlib.Path('reference_Repair.pdb').write_text(data)
        pathlib.Path('reference_Repair.fxout').write_text('Synthetic repair log\n')
else:
    assert args['numberOfRuns'] == '1' and args['out-pdb'] == 'true'
    lines = pathlib.Path(args['mutant-file']).read_text().splitlines()
    header = 'Pdb\ttotal energy\tBackbone Hbond\n'
    rows = []
    three = {'V':'VAL','S':'SER'}
    for i, mutation in enumerate(lines, 1):
        if mode == 'truncated' and i == len(lines): continue
        name = f'reference_Repair_{i}.pdb'
        rows.append(f'{name}\t{i}.25\t-0.1\n')
        source = pathlib.Path(args['pdb']).read_text()
        mutant = []
        tokens = mutation.rstrip(';').split(',')
        for line in source.splitlines(keepends=True):
            if line.startswith('ATOM') and mode != 'unmutated':
                for token in tokens:
                    if line[21] == token[1] and int(line[22:26]) == int(token[2:-1]):
                        line = line[:17]+three[token[-1]]+line[20:]
            mutant.append(line)
        pathlib.Path(name).write_text(''.join(mutant))
        shutil.copyfile(args['pdb'], 'WT_' + name)
    if mode == 'wrong_mapping': rows.reverse()
    pathlib.Path('Dif_reference_Repair.fxout').write_text(header + ''.join(rows))
    for prefix in ('Raw_', 'Average_', 'PdbList_'):
        if mode != 'missing_raw' or prefix != 'Raw_':
            pathlib.Path(prefix + 'reference_Repair.fxout').write_text('SYNTHETIC fixture\n')
'''


def atom(serial, name, residue, number):
    return f"ATOM  {serial:5d} {name:>4} {residue:>3} A{number:4d}    {1.0:8.3f}{2.0:8.3f}{3.0:8.3f}{1.0:6.2f}{20.0:6.2f}          {name[0]:>2}\n"


class FoldXTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.fasta = self.root / 'reference.fasta'
        self.fasta.write_text('>test\nAC\n')
        self.pdb = self.root / 'reference.pdb'
        self.pdb.write_text(''.join(atom((res-1)*4+i,name,aa,res) for res,aa in ((1,'ALA'),(2,'CYS')) for i,name in enumerate(('N','CA','C','O'),1))+'END\n')
        self.config = {
            'schema_version':1,'study_id':'test-study','comparison_type':'point_mutation','synthetic':True,
            'reference':{'fasta':str(self.fasta),'pdb':str(self.pdb),'chain_map':{'A':'test'}},
            'variants':[{'id':'A1V','mutations':[{'chain':'A','residue_number':1,'from':'A','to':'V'}]}, {'id':'C2S','mutations':[{'chain':'A','residue_number':2,'from':'C','to':'S'}]}],
            'foldx':{'runs':2,'version':'SYNTHETIC-FAKE-NOT-FOLDX','timeout_seconds':10}}

    def tearDown(self): self.temp.cleanup()

    def fake(self, mode='ok'):
        path = self.root / ('fake foldx ; ' + mode)
        path.write_text('#!'+sys.executable+'\n'+FAKE.replace('MODE',repr(mode)),encoding='utf-8')
        path.chmod(0o700)
        return path

    def test_plan_snapshots_mapping_original_unchanged_and_exclusive_outdir(self):
        before = self.pdb.read_bytes()
        original = copy.deepcopy(self.config)
        plan = plan_study(self.config,self.root/'plan')
        self.assertEqual(self.config,original)
        self.assertEqual(self.pdb.read_bytes(),before)
        self.assertEqual(plan['variants'][0]['mutation_line'],'AA1V;')
        self.assertEqual(plan['reference']['pdb_sha256'],hashlib.sha256(before).hexdigest())
        self.assertEqual((self.root/'plan/individual_list.txt').read_text(),'AA1V;\nCA2S;\n')
        self.assertTrue((self.root/'plan/prepared/reference.fasta').is_file())
        self.assertTrue(plan['synthetic'])
        self.assertIn('--numberOfRuns=1',plan['commands']['build_per_repeat'])
        with self.assertRaisesRegex(ValueError,'already exists'): plan_study(self.config,self.root/'plan')

    def test_unsupported_comparisons_no_output(self):
        for comparison in ('linker','domain_order','cross_protein',[],{}):
            self.config['comparison_type']=comparison
            out=self.root/'unsupported'
            with self.subTest(comparison=comparison),self.assertRaises(ValueError): plan_study(self.config,out)
            self.assertFalse(out.exists())

    def test_strict_types_ranges_unknown_fields(self):
        for field,value in (('runs',True),('runs',0),('runs',1.5),('pH',float('nan')),('pH',15),('temperature',0),('ion_strength',-.1),('timeout_seconds',True)):
            config=copy.deepcopy(self.config);config['foldx'][field]=value
            with self.subTest(field=field,value=value),self.assertRaises(ValueError): validate_config(config)
        self.config['foldx']['ignored_flag']=1
        with self.assertRaisesRegex(ValueError,'unknown fields'): validate_config(self.config)

    def test_mutations_and_wt(self):
        for field,value in (('residue_number',True),('residue_number',-1),('from','Z'),('to','A'),('chain','B')):
            config=copy.deepcopy(self.config);config['variants'][0]['mutations'][0][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError): validate_config(config)
        self.config['variants'][0]['mutations'][0]['from']='G'
        with self.assertRaisesRegex(ValueError,'WT mismatch'): plan_study(self.config,self.root/'mismatch')

    def test_license_after_offline_plan(self):
        out=self.root/'unlicensed'
        with self.assertRaisesRegex(ValueError,'license_confirmed=True'): run_study(self.config,out,self.fake())
        self.assertTrue((out/'plan.json').is_file())
        self.assertEqual(json.loads((out/'run_manifest.json').read_text())['status'],'failed')
        self.assertFalse((out/'repair').exists())

    def test_synthetic_success_full_provenance_and_collection(self):
        out=self.root/'success';binary=self.fake()
        manifest=run_study(self.config,out,binary,license_confirmed=True)
        self.assertEqual(manifest['status'],'completed')
        self.assertEqual(len(manifest['runs']),2)
        self.assertEqual(len(manifest['steps']),3)
        self.assertEqual(manifest['executable']['sha256'],hashlib.sha256(binary.read_bytes()).hexdigest())
        for run in manifest['runs']:
            self.assertEqual(run['diff_sha256'],hashlib.sha256((out/run['diff_file']).read_bytes()).hexdigest())
        self.assertIn('SYNTHETIC',(out/'repair/stdout.log').read_text())
        self.assertIn('synthetic stderr',(out/'repair/stderr.log').read_text())
        self.assertTrue(collect_results(out)['synthetic'])

    def test_real_manifest_contract_using_fake_process_is_archive_portable(self):
        # Exercise strict non-synthetic provenance validation only. The process,
        # coordinates and numbers remain fake, temporary software-test data.
        self.config['synthetic'] = False
        out = self.root / 'original-execution'
        run_study(self.config, out, self.fake(), license_confirmed=True)
        archive = self.root / 'moved-archive'
        out.rename(archive)
        self.fasta.unlink()
        self.pdb.unlink()
        result = collect_results(archive)
        self.assertFalse(result['synthetic'])
        self.assertEqual(result['study_id'], 'test-study')

    def test_external_errors_incomplete_outputs_identity_mapping(self):
        for mode in ('exit','missing_repaired','changed_repair','missing_raw','truncated','wrong_mapping','unmutated'):
            out=self.root/mode
            with self.subTest(mode=mode),self.assertRaises(ValueError): run_study(self.config,out,self.fake(mode),license_confirmed=True)
            manifest=json.loads((out/'run_manifest.json').read_text())
            self.assertEqual(manifest['status'],'failed')
            self.assertIn('error',manifest)
            self.assertTrue((out/'repair/stdout.log').is_file())

    def test_timeout_and_missing_executable(self):
        self.config['foldx']['timeout_seconds']=1
        out=self.root/'timeout-run'
        with self.assertRaisesRegex(ValueError,'exceeded 1 seconds'): run_study(self.config,out,self.fake('timeout'),license_confirmed=True)
        self.assertTrue(json.loads((out/'run_manifest.json').read_text())['steps'][0]['timed_out'])
        with self.assertRaisesRegex(ValueError,'executable'): run_study(self.config,self.root/'missing',self.root/'absent',license_confirmed=True)

    def test_resources_isolated_hashed_collision_blocked(self):
        resources=self.root/'licensed-resources';resources.mkdir()
        original=resources/'rotabase.txt';original.write_text('SYNTHETIC resource\n')
        self.config['foldx']['resources_dir']=str(resources)
        out=self.root/'resources'
        manifest=run_study(self.config,out,self.fake(),license_confirmed=True)
        self.assertEqual(manifest['resources'][0]['sha256'],hashlib.sha256(original.read_bytes()).hexdigest())
        for name in ('repair','replicates/run_001','replicates/run_002'): self.assertEqual((out/name/'rotabase.txt').read_bytes(),original.read_bytes())
        (resources/'reference.pdb').write_text('collision')
        with self.assertRaisesRegex(ValueError,'collide'): plan_study(self.config,self.root/'collision')

    def test_multimutant_order_and_duplicate_design(self):
        self.config['variants']=[{'id':'double','mutations':[self.config['variants'][1]['mutations'][0],self.config['variants'][0]['mutations'][0]]}]
        plan=plan_study(self.config,self.root/'double')
        self.assertEqual(plan['variants'][0]['mutation_line'],'AA1V,CA2S;')
        self.config['variants'].append({'id':'same','mutations':list(reversed(self.config['variants'][0]['mutations']))})
        with self.assertRaisesRegex(ValueError,'duplicate mutation design'): validate_config(self.config)

if __name__=='__main__': unittest.main()
