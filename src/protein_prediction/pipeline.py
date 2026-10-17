"""Dependency-free, auditable FoldX point-substitution workflow."""
import argparse
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import sys

AA = 'ACDEFGHIKLMNPQRSTVWY'
THREE = dict(zip('ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(), AA))

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def fasta(path):
    records = {}
    name = None
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            name = line[1:].split()[0] if line[1:].strip() else ''
            if not name or name in records:
                raise ValueError('Empty or duplicate FASTA identifier')
            records[name] = ''
        else:
            if name is None:
                raise ValueError('Sequence before FASTA header')
            records[name] += line.upper()
    if not records or any(not s or set(s) - set(AA) for s in records.values()):
        raise ValueError('FASTA must contain nonempty sequences using the 20 canonical amino acids')
    return records

def pdb(path):
    """Strict preparation gate. Preserve numbering/coordinates; never invent missing atoms."""
    residues = {}
    atoms = set()
    models = 0
    for line in Path(path).read_text().splitlines():
        if line.startswith('MODEL '):
            models += 1
        if line.startswith('HETATM'):
            raise ValueError('Resolve HETATM/cofactors/waters explicitly before preparation')
        if not line.startswith('ATOM  '):
            continue
        if len(line) < 54 or line[16] != ' ' or line[26] != ' ':
            raise ValueError('Resolve alternate locations/insertion codes before preparation')
        chain, resname, atom = line[21], line[17:20], line[12:16].strip()
        if not chain.isalnum() or resname not in THREE:
            raise ValueError('Use alphanumeric chain IDs and canonical residues')
        number = int(line[22:26])
        xyz = [float(line[i:i+8]) for i in (30, 38, 46)]
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError('Nonfinite PDB coordinates')
        key = (chain, number)
        if (key, atom) in atoms:
            raise ValueError('Duplicate atom')
        atoms.add((key, atom))
        if key in residues and residues[key]['aa'] != THREE[resname]:
            raise ValueError('Conflicting residue identities')
        residues.setdefault(key, {'aa': THREE[resname], 'atoms': set()})['atoms'].add(atom)
    if models > 1 or not residues:
        raise ValueError('Provide one nonempty PDB model')
    for key, value in residues.items():
        if not {'N', 'CA', 'C', 'O'} <= value['atoms']:
            raise ValueError(f'Missing backbone atoms at {key}')
    return residues

def prepare(config_path, output, execute=False, binary=None):
    config_path = Path(config_path).resolve()
    cfg = json.loads(config_path.read_text())
    if cfg.get('comparison') != 'point_mutation':
        raise ValueError('BuildModel supports this workflow only for same-backbone point substitutions; see docs/comparisons.md')
    count = cfg.get('repeats', 5)
    if type(count) is not int or count < 2 or count > 100:
        raise ValueError('repeats must be an integer from 2 to 100')
    cfg['repeats'] = count
    for k, lo, hi in [('pH', 0, 14), ('temperature', 200, 400), ('ionStrength', 0, 5)]:
        val = cfg.get(k, {'pH': 7, 'temperature': 298, 'ionStrength': .05}[k])
        if not isinstance(val, (int, float)) or not math.isfinite(val) or not lo <= val <= hi:
            raise ValueError(f'Invalid {k}')
    source = (config_path.parent / cfg['pdb']).resolve()
    seqpath = (config_path.parent / cfg['fasta']).resolve()
    residues, sequences = pdb(source), fasta(seqpath)
    mapping = cfg['chains']
    chains = {c for c, n in residues}
    if set(mapping) != chains:
        raise ValueError('chains mapping must cover every PDB chain exactly')
    for chain, record in mapping.items():
        observed = ''.join(v['aa'] for (c, n), v in residues.items() if c == chain)
        if sequences.get(record) != observed:
            raise ValueError(f'FASTA/PDB sequence mismatch for chain {chain}; inspect unresolved residues and numbering')
    variants = cfg.get('variants', [])
    if not variants:
        raise ValueError('At least one variant required')
    names = set()
    mutation_lines = []
    for variant in variants:
        name = variant['id']
        if not re.fullmatch(r'[A-Za-z0-9_-]+', name) or name in names:
            raise ValueError('Variant IDs must be unique safe identifiers')
        names.add(name)
        seen, codes = set(), []
        for mutation in variant['mutations']:
            match = re.fullmatch(r'(['+AA+r'])([A-Za-z0-9])(-?\d+)(['+AA+r'])', mutation)
            if not match:
                raise ValueError(f'Invalid mutation: {mutation}')
            wt, chain, number, alt = match.groups()
            key = chain, int(number)
            if key in seen or key not in residues or residues[key]['aa'] != wt or wt == alt:
                raise ValueError(f'Duplicate, missing, WT mismatch or no-op mutation: {mutation}')
            seen.add(key)
            codes.append(mutation)
        if not codes:
            raise ValueError('Empty mutation list')
        mutation_lines.append(','.join(codes)+';')
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(source, out/'input.pdb')
    shutil.copyfile(seqpath, out/'input.fasta')
    (out/'individual_list.txt').write_text('\n'.join(mutation_lines)+'\n')
    params = [f'--{k}={cfg.get(k, default)}' for k, default in [('pH', 7), ('temperature', 298), ('ionStrength', .05)]]
    commands = [ ['--command=RepairPDB', '--pdb=input.pdb']+params,
                 ['--command=BuildModel', '--pdb=input_Repair.pdb', '--mutant-file=individual_list.txt', f'--numberOfRuns={count}']+params ]
    manifest = {'schema': 1, 'status': 'planned', 'config': cfg, 'commands': commands,
                'inputs': {'pdb_sha256': digest(source), 'fasta_sha256': digest(seqpath)},
                'python': sys.version, 'platform': platform.platform(),
                'notice': 'Repeated runs measure algorithmic spread, not biological uncertainty.'}
    def save():
        (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    save()
    if execute:
        try:
            exe = shutil.which(binary or os.environ.get('FOLDX_BIN', 'foldx'))
            if not exe:
                raise ValueError('FoldX not found; supply licensed executable via --foldx or FOLDX_BIN')
            exe = str(Path(exe).resolve())
            manifest['binary_sha256'] = digest(exe)
            manifest['binary'] = exe
            help_result = subprocess.run([exe, '--help'], capture_output=True, text=True, timeout=60)
            (out/'foldx-help.txt').write_text(help_result.stdout+help_result.stderr)
            # User-supplied ancillary files remain private and are checksummed.
            ancillary = cfg.get('rotabase')
            if ancillary:
                rotabase = (config_path.parent / ancillary).resolve()
                shutil.copyfile(rotabase, out/'rotabase.txt')
                manifest['rotabase_sha256'] = digest(rotabase)
            for i, command in enumerate(commands):
                with (out/f'step-{i+1}.log').open('w') as log:
                    subprocess.run([exe]+command, cwd=out, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=cfg.get('timeout_seconds', 3600))
                if i == 0:
                    repaired = pdb(out/'input_Repair.pdb')
                    if {k:v['aa'] for k,v in repaired.items()} != {k:v['aa'] for k,v in residues.items()}:
                        raise ValueError('RepairPDB changed residue mapping')
            files = list(out.glob('Dif_*.fxout'))
            if len(files) != 1:
                raise ValueError('Expected exactly one Dif_*.fxout')
            summarize(files[0], out/'manifest.json', out/'report', config=cfg)
            manifest['status'] = 'completed'
            manifest['outputs'] = {p.name:digest(p) for p in out.iterdir() if p.is_file() and p.name != 'manifest.json'}
        except Exception as error:
            manifest['status'] = 'failed'
            manifest['error'] = str(error)
            raise
        finally:
            save()
    return manifest

def parse_dif(path):
    if not Path(path).name.startswith('Dif_'):
        raise ValueError('Only per-run Dif_*.fxout accepted; Average/Raw are not paired per-run ddG')
    header = None
    rows = []
    for line in Path(path).read_text().splitlines():
        cells = [x.strip() for x in line.split('\t')]
        if cells[0].lower() == 'pdb' and 'total energy' in [x.lower() for x in cells]:
            header = [x.lower() for x in cells]
            continue
        if header is None or not line.strip():
            continue
        if len(cells) != len(header):
            raise ValueError('Malformed fxout row')
        value = float(cells[header.index('total energy')])
        if not math.isfinite(value):
            raise ValueError('Nonfinite energy')
        rows.append((cells[0], value))
    if not rows:
        raise ValueError('Missing header or no result rows')
    return rows

def summarize(path, manifest_path, output, config=None):
    manifest = json.loads(Path(manifest_path).read_text())
    cfg = config or manifest['config']
    if cfg['comparison'] != 'point_mutation':
        raise ValueError('Cannot rank different architectures using point-mutation ddG')
    variants, repeats = cfg['variants'], cfg['repeats']
    groups = {i:[] for i in range(1, len(variants)+1)}
    ids = set()
    # FoldX BuildModel names: <pdb-stem>_<variant-index>_<run-index>.pdb
    for label, value in parse_dif(path):
        match = re.fullmatch(r'.+_(\d+)_(\d+)(?:\.pdb)?', label)
        if not match:
            raise ValueError(f'Unsupported FoldX row name {label}; do not infer variant identity from row order')
        variant, run = map(int, match.groups())
        if variant not in groups or not 1 <= run <= repeats or (variant, run) in ids:
            raise ValueError('Unexpected or duplicate variant/run')
        ids.add((variant, run))
        groups[variant].append(value)
    if any(len(values) != repeats for values in groups.values()):
        raise ValueError('Incomplete replicate coverage')
    result = []
    for i, variant in enumerate(variants, 1):
        values = groups[i]
        result.append({'variant':variant['id'], 'mutations':','.join(variant['mutations']), 'n':len(values),
                       'mean_ddg':statistics.mean(values), 'sd':statistics.stdev(values), 'min':min(values), 'max':max(values)})
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    (out/'summary.json').write_text(json.dumps({'data_kind':cfg.get('data_kind','user_supplied'), 'source_sha256':digest(path), 'units':'kcal/mol', 'results':result}, indent=2)+'\n')
    with (out/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)
    trs = ''.join('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in row.values())+'</tr>' for row in result)
    scale = max(1, max(abs(r['min']) for r in result), max(abs(r['max']) for r in result))
    shapes = ''
    for i, row in enumerate(result):
        y = 35+i*36
        x = 380+row['mean_ddg']/scale*180
        lo, hi = 380+row['min']/scale*180, 380+row['max']/scale*180
        shapes += f'<text x="10" y="{y+4}">{html.escape(row["variant"])}</text><line x1="{lo}" x2="{hi}" y1="{y}" y2="{y}" stroke="#234d63"/><circle cx="{x}" cy="{y}" r="4" fill="#234d63"/>'
    height = 70+36*len(result)
    page = f'''<!doctype html><html lang="zh"><meta charset="utf-8"><title>Stability analysis</title>
<style>body{{background:#faf8f0;color:#222;font:16px system-ui;max-width:1100px;margin:40px auto;padding:20px}}table{{border-collapse:collapse;width:100%}}td,th{{padding:10px;border-bottom:1px solid #ccc;text-align:left}}svg{{width:100%}}</style>
<h1>抗体工程化蛋白稳定性分析</h1><p>数据类型：{html.escape(cfg.get('data_kind','user_supplied'))}。同一骨架点突变；ΔΔG = G(mutant) − G(WT)，正值表示预测去稳定化，负值表示预测稳定化，单位 kcal/mol。</p>
<p>点为均值，横线为运行最小值至最大值；表中 SD 为算法重复的样本标准差。重复计算不能代表模型误差，也不能证明热稳定性或表达表现。</p>
<svg role="img" aria-label="Mean ddG and run range" viewBox="0 0 600 {height}"><line x1="380" x2="380" y1="10" y2="{height-20}" stroke="#aaa"/>{shapes}<text x="200" y="{height-4}">−{scale:.2f}</text><text x="375" y="{height-4}">0</text><text x="550" y="{height-4}">{scale:.2f}</text></svg>
<table><tr>{''.join('<th>'+html.escape(k)+'</th>' for k in result[0])}</tr>{trs}</table>
<h2>结果判读</h2><p>结合结构来源、突变位点置信度和多结构一致性评估。linker、VH/VL 顺序、拓扑或蛋白长度变化须采用独立结构评估，不以本图排序。实验验证建议 nanoDSF/DSC、SEC 单体比例及应力后聚集。</p><p>结果文件 SHA256：{digest(path)}</p></html>'''
    (out/'index.html').write_text(page)
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('validate-fasta'); p.add_argument('path')
    p = sub.add_parser('inspect-pdb'); p.add_argument('path')
    for command in ('plan','run'):
        p = sub.add_parser(command); p.add_argument('config'); p.add_argument('--out',required=True); p.add_argument('--foldx')
    p = sub.add_parser('report'); p.add_argument('dif'); p.add_argument('--manifest',required=True); p.add_argument('--out',required=True)
    sub.add_parser('doctor')
    args = parser.parse_args()
    try:
        if args.command == 'validate-fasta':
            print(json.dumps({k:{'length':len(v), 'cysteines':v.count('C')} for k,v in fasta(args.path).items()},indent=2))
        elif args.command == 'inspect-pdb':
            data=pdb(args.path); print(json.dumps({'residues':len(data),'chains':sorted({c for c,n in data})}))
        elif args.command in ('plan','run'):
            print(json.dumps(prepare(args.config,args.out,args.command=='run',args.foldx),indent=2))
        elif args.command == 'report':
            print(json.dumps(summarize(args.dif,args.manifest,args.out),indent=2))
        else:
            print(json.dumps({'python':sys.version,'platform':platform.platform(),'foldx':shutil.which(os.environ.get('FOLDX_BIN','foldx')),'license':'User must verify entitlement; executable presence does not establish license'},indent=2))
    except (ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        parser.exit(2,f'Error: {error}\n')

if __name__ == '__main__':
    main()
