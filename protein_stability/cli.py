"""Offline-first command line interface. No service credentials are required."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import sys

from . import __version__


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key in study configuration: {key}")
        result[key] = value
    return result


def load_config(path):
    """Resolve file paths relative to the manifest, never the working directory."""
    path = Path(path).resolve()
    config = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_unique_object)
    if not isinstance(config, dict):
        raise ValueError("Study config must be a JSON object")
    reference = config.get("reference")
    if isinstance(reference, dict):
        for key in ("fasta", "pdb"):
            if isinstance(reference.get(key), str):
                reference[key] = str((path.parent / reference[key]).resolve())
    foldx = config.get("foldx")
    if isinstance(foldx, dict) and isinstance(foldx.get("resources_dir"), str):
        foldx["resources_dir"] = str((path.parent / foldx["resources_dir"]).resolve())
    return config


def doctor(executable=None):
    binary = executable or os.environ.get("FOLDX_BINARY") or "foldx"
    found = shutil.which(binary)
    return {
        "schema_version": 1,
        "package_version": __version__,
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "architecture": platform.machine(),
        "runtime_dependencies": "Python standard library only",
        "offline_modules_ready": sys.version_info >= (3, 10),
        "foldx_executable": found,
        "foldx_detected": found is not None,
        "foldx_execution_verified": False,
        "foldx_license_verified": False,
        "container_engine": shutil.which("docker") or shutil.which("podman"),
        "optional_modules": {name: importlib.util.find_spec(name) is not None
                             for name in ("numpy", "matplotlib", "Bio")},
        "next_step": "Supply a licensed platform-compatible FoldX binary and resources, "
                     "complete structures and exact chain FASTA mapping; then run a local smoke test.",
    }


def emit(data, output=None):
    text = json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    if output:
        destination = Path(output)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(text, encoding="utf-8")
    print(text, end="")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("doctor", help="Check environment without executing FoldX")
    p.add_argument("--foldx")
    p.add_argument("--output")
    p = commands.add_parser("validate-fasta", help="Validate canonical amino-acid FASTA")
    p.add_argument("fasta")
    p.add_argument("--output")
    p = commands.add_parser("inspect-pdb", help="Inspect PDB sequence, numbering and atom completeness")
    p.add_argument("pdb")
    p.add_argument("--output")
    p = commands.add_parser("prepare-pdb", help="Extract explicitly mapped chains and verify sequence")
    p.add_argument("pdb")
    p.add_argument("--fasta", required=True)
    p.add_argument("--chain-map", required=True, help='JSON object, e.g. \'{"A":"reference"}\'')
    p.add_argument("--output", required=True, help="Prepared PDB destination")
    p.add_argument("--audit", help="Optional JSON audit destination")
    for command in ("plan", "run"):
        p = commands.add_parser(command, help="Plan or execute a point-mutation study")
        p.add_argument("config")
        p.add_argument("--out", required=True, help="New output directory (must not exist)")
        if command == "run":
            p.add_argument("--foldx", default=os.environ.get("FOLDX_BINARY", "foldx"))
            p.add_argument("--license-confirmed", action="store_true",
                           help="Confirm authorization to use your supplied FoldX installation")
    p = commands.add_parser("report", help="Parse a completed run and write an HTML/CSV/JSON report")
    p.add_argument("run_dir")
    p.add_argument("--out", required=True)
    p = commands.add_parser("demo", help="Generate synthetic data and test reporting without FoldX")
    p.add_argument("--out", required=True)
    p = commands.add_parser("comparison-template", help="Create an unscored structural-comparison report")
    p.add_argument("--type", choices=("linker", "domain_order", "cross_protein"), required=True)
    p.add_argument("--study-id", default="structural-comparison")
    p.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    try:
        from .inputs import inspect_fasta, inspect_pdb, prepare_pdb, read_fasta
        from .foldx import plan_study, run_study
        from .results import collect_results
        from .report import make_demo, write_report
        if args.command == "doctor":
            emit(doctor(args.foldx), args.output)
        elif args.command == "validate-fasta":
            emit(inspect_fasta(args.fasta), args.output)
        elif args.command == "inspect-pdb":
            emit(inspect_pdb(args.pdb), args.output)
        elif args.command == "prepare-pdb":
            mapping = json.loads(args.chain_map)
            emit(prepare_pdb(args.pdb, args.output, mapping, read_fasta(args.fasta)), args.audit)
        elif args.command in ("plan", "run"):
            config = load_config(args.config)
            if args.command == "plan":
                emit(plan_study(config, args.out))
            else:
                emit(run_study(config, args.out, args.foldx,
                               license_confirmed=args.license_confirmed))
        elif args.command == "report":
            emit(write_report(collect_results(args.run_dir), args.out))
        elif args.command == "demo":
            emit(make_demo(args.out))
        elif args.command == "comparison-template":
            emit(write_report({
                "schema_version": 1, "study_id": args.study_id,
                "comparison_type": args.type, "synthetic": False,
                "status": "requires_structural_ensemble", "variants": [],
                "provenance": {"template": True},
                "warnings": ["尚未提供模型与实验数据；本模板不产生 ΔΔG 或稳定性排名。"],
            }, args.out))
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0
