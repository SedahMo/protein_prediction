"""Strict parsing and scoped aggregation of FoldX mutation differences.

Only Dif_/Diff_ files produced by BuildModel describe mutation differences.
Raw_ and Average_ files contain other quantities and are deliberately rejected.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
from typing import Any


class ResultsError(ValueError):
    """An output or manifest is incomplete, malformed, or out of scope."""


def parse_diff(path: str | Path) -> list[dict[str, Any]]:
    """Parse a complete tab-delimited BuildModel Dif_/Diff_ output.

    Preamble lines before the ``Pdb\ttotal energy`` header are permitted. Every
    subsequent nonblank line must be a full row with finite numeric terms.
    ``total energy`` is returned as mutant-minus-WT ΔΔG in kcal/mol; a positive
    value means destabilization in the FoldX model.
    """
    path = Path(path)
    if not path.name.startswith(("Dif_", "Diff_")) or path.suffix != ".fxout":
        raise ResultsError("ΔΔG parsing requires a Dif_*.fxout or Diff_*.fxout file; Raw_/Average_ are not mutation differences")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            records = list(csv.reader(handle, delimiter="\t", strict=True))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ResultsError(f"Cannot read FoldX difference file {path}: {exc}") from exc
    header: list[str] | None = None
    terminal_tab = False
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for line_no, cells in enumerate(records, 1):
        if not cells or all(not cell.strip() for cell in cells):
            continue
        cells = [cell.strip() for cell in cells]
        if header is None:
            if cells[0] != "Pdb":
                continue
            if len(cells) < 2 or cells[1] != "total energy":
                raise ResultsError(f"{path.name}:{line_no}: expected tab-delimited Pdb and total energy header")
            # Some FoldX versions terminate every table line with one tab.
            # Accept that convention consistently, without hiding empty terms.
            terminal_tab = cells[-1] == ""
            if terminal_tab:
                cells = cells[:-1]
            normalized = [" ".join(cell.casefold().split()) for cell in cells]
            if any(not cell for cell in cells) or len(normalized) != len(set(normalized)):
                raise ResultsError(f"{path.name}:{line_no}: empty or duplicate column header")
            header = cells
            continue
        if terminal_tab:
            if cells[-1] != "":
                raise ResultsError(f"{path.name}:{line_no}: inconsistent terminal-tab convention")
            cells = cells[:-1]
        if len(cells) != len(header) or any(not cell for cell in cells):
            raise ResultsError(f"{path.name}:{line_no}: incomplete or malformed row: expected {len(header)} nonempty tab-separated cells, got {len(cells)}")
        model_id = cells[0]
        if model_id == "Pdb" or model_id in seen:
            raise ResultsError(f"{path.name}:{line_no}: duplicate header or model id {model_id!r}")
        values = []
        for column, value in zip(header[1:], cells[1:]):
            try:
                number = float(value)
            except ValueError as exc:
                raise ResultsError(f"{path.name}:{line_no}: nonnumeric {column!r}: {value!r}") from exc
            if not math.isfinite(number):
                raise ResultsError(f"{path.name}:{line_no}: nonfinite {column!r}")
            if re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", value) is None:
                raise ResultsError(f"{path.name}:{line_no}: malformed numeric {column!r}: {value!r}")
            values.append(number)
        seen.add(model_id)
        rows.append({"model_id": model_id, "ddg_kcal_mol": values[0], "terms": dict(zip(header[2:], values[1:]))})
    if header is None:
        raise ResultsError(f"{path.name}: no Pdb\ttotal energy header found")
    if not rows:
        raise ResultsError(f"{path.name}: header found but no data rows; output may be truncated")
    return rows


def _unique_json(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ResultsError(f"Duplicate JSON key in run manifest: {key}")
        result[key] = value
    return result


def _positive_int(value: Any, label: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ResultsError(f"{label} must be a positive integer")
    return value


def _text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ResultsError(f"{label} must be a nonempty string")
    return value


def _relative_path(root: Path, value: Any, label: str) -> Path:
    value = _text(value, label)
    if Path(value).is_absolute():
        raise ResultsError(f"{label} must be relative to the manifest directory")
    path = (root / value).resolve()
    if not path.is_relative_to(root):
        raise ResultsError(f"{label} escapes the manifest directory")
    return path


def _model_key(value: str) -> str:
    """FoldX versions differ in whether model labels retain '.pdb'."""
    return value[:-4] if value.endswith(".pdb") else value


def _check_artifact(root: Path, relative: str, expected: Any, label: str) -> None:
    if not isinstance(expected, str) or len(expected) != 64 or any(char not in "0123456789abcdef" for char in expected):
        raise ResultsError(f"{label} must record a lowercase SHA-256 digest")
    path = _relative_path(root, relative, label)
    try:
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        raise ResultsError(f"Complete run folder required: cannot read {relative}: {exc}") from exc
    if actual != expected:
        raise ResultsError(f"{label}: SHA-256 mismatch for {relative}")


def _verify_genuine_run(root: Path, manifest: dict[str, Any]) -> None:
    """Verify local snapshots and the immutable plan, never original input paths."""
    reference = manifest["reference"]
    for relative, key in (("prepared/source.pdb", "pdb_sha256"), ("prepared/reference.fasta", "fasta_sha256"), ("prepared/reference.pdb", "prepared_pdb_sha256"), ("repair/reference_Repair.pdb", "repaired_pdb_sha256")):
        _check_artifact(root, relative, reference.get(key), f"reference.{key}")
    runtime = manifest.get("runtime")
    if not isinstance(runtime, dict):
        raise ResultsError("Genuine run requires runtime.plan_sha256 and the complete local plan")
    _check_artifact(root, "plan.json", runtime.get("plan_sha256"), "runtime.plan_sha256")
    try:
        plan = json.loads((root / "plan.json").read_text(encoding="utf-8-sig"), object_pairs_hook=_unique_json)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ResultsError(f"Cannot read verified plan.json: {exc}") from exc
    if not isinstance(plan, dict) or type(plan.get("schema_version")) is not int or plan.get("schema_version") != 1:
        raise ResultsError("Genuine run plan has an invalid schema")
    for key in ("study_id", "comparison_type", "synthetic", "variants"):
        if plan.get(key) != manifest.get(key):
            raise ResultsError(f"Plan and run manifest disagree on {key}")
    if plan.get("reference") != {key: value for key, value in reference.items() if key != "repaired_pdb_sha256"}:
        raise ResultsError("Plan and run manifest disagree on the reference backbone/sequence")
    config = plan.get("config")
    if not isinstance(config, dict) or not isinstance(config.get("foldx"), dict):
        raise ResultsError("Genuine plan must retain the FoldX configuration")
    if any(manifest["foldx"].get(key) != value for key, value in config["foldx"].items()):
        raise ResultsError("Plan and run manifest disagree on FoldX settings")
    _check_artifact(root, "individual_list.txt", plan.get("mutation_file_sha256"), "plan.mutation_file_sha256")


def _verify_step_artifacts(root: Path, steps: Any, runs: list[dict[str, Any]]) -> None:
    """Check recorded outputs while allowing a run folder to be moved/archived."""
    if steps is None:
        return
    if not isinstance(steps, list):
        raise ResultsError("Recorded steps must be a list")
    directories = {run["replicate"]: run["directory"] for run in runs}
    for step in steps:
        if not isinstance(step, dict) or step.get("returncode") != 0:
            raise ResultsError("A completed run contains an unsuccessful or malformed execution step")
        if step.get("name") == "RepairPDB":
            cwd = root / "repair"
        elif step.get("name") == "BuildModel" and type(step.get("replicate")) is int and step["replicate"] in directories:
            cwd = _relative_path(root, directories[step["replicate"]], "run.directory")
        else:
            raise ResultsError("Unknown or unscoped recorded execution step")
        files = step.get("files")
        if not isinstance(files, list):
            raise ResultsError("Recorded execution step is missing output artifacts")
        for artifact in files:
            if not isinstance(artifact, dict):
                raise ResultsError("Recorded artifacts must be objects")
            path = _relative_path(cwd.resolve(), artifact.get("path"), "step.artifact.path")
            _check_artifact(root, str(path.relative_to(root)), artifact.get("sha256"), "step.artifact.sha256")


def collect_results(run_dir: str | Path) -> dict[str, Any]:
    """Collect exactly one completed point-mutation study and reference.

    Missing repeats, unknown models, inconsistent mutation mappings, duplicate
    outputs, and mismatched study/reference scopes fail instead of generating a
    partial ranking. FoldX row order is checked against ``foldx_index`` as well
    as each expected model label.
    """
    root = Path(run_dir).resolve()
    manifest_path = root / "run_manifest.json"
    try:
        manifest_bytes = manifest_path.read_bytes()
        manifest = json.loads(manifest_bytes.decode("utf-8-sig"), object_pairs_hook=_unique_json)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ResultsError(f"Cannot read {manifest_path}: {exc}") from exc
    if not isinstance(manifest, dict) or type(manifest.get("schema_version")) is not int or manifest.get("schema_version") != 1:
        raise ResultsError("Unsupported run manifest; schema_version must be 1")
    study_id = _text(manifest.get("study_id"), "study_id")
    if manifest.get("comparison_type") != "point_mutation":
        raise ResultsError("FoldX ΔΔG aggregation is restricted to point_mutation studies on one reference backbone")
    if manifest.get("status") not in ("complete", "completed"):
        raise ResultsError("Run is incomplete: manifest status must be completed before aggregation")
    if not isinstance(manifest.get("synthetic"), bool):
        raise ResultsError("Manifest must explicitly declare synthetic: true or false")
    reference = manifest.get("reference")
    if not isinstance(reference, dict):
        raise ResultsError("Manifest must identify one reference structure and sequence")
    for key in ("pdb_sha256", "fasta_sha256", "prepared_pdb_sha256"):
        digest = reference.get(key)
        if not isinstance(digest, str) or len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise ResultsError(f"reference.{key} must be a lowercase SHA-256 digest")
    foldx = manifest.get("foldx")
    if not isinstance(foldx, dict):
        raise ResultsError("Manifest is missing FoldX settings")
    expected_n = _positive_int(foldx.get("runs"), "foldx.runs")
    if not manifest["synthetic"]:
        _verify_genuine_run(root, manifest)
    variants = manifest.get("variants")
    if not isinstance(variants, list) or not variants:
        raise ResultsError("Manifest must contain at least one variant")
    variant_ids: list[str] = []
    variant_by_index: dict[int, dict[str, Any]] = {}
    for variant in variants:
        if not isinstance(variant, dict):
            raise ResultsError("Every variant must be an object")
        variant_id = _text(variant.get("id"), "variant.id")
        index = _positive_int(variant.get("foldx_index"), "variant.foldx_index")
        _text(variant.get("mutation_line"), "variant.mutation_line")
        if variant_id in variant_ids or index in variant_by_index:
            raise ResultsError("Duplicate variant id or FoldX index")
        variant_ids.append(variant_id)
        variant_by_index[index] = variant
    if set(variant_by_index) != set(range(1, len(variants) + 1)):
        raise ResultsError("Variant FoldX indices must be consecutive from 1")
    runs = manifest.get("runs")
    if not isinstance(runs, list) or len(runs) != expected_n:
        raise ResultsError(f"Expected {expected_n} replicate runs; every replicate must be present")
    measurements: dict[str, list[dict[str, Any]]] = {variant_id: [] for variant_id in variant_ids}
    seen_replicates: set[int] = set()
    seen_files: set[Path] = set()
    provenance_files: list[dict[str, Any]] = []
    for run in runs:
        if not isinstance(run, dict):
            raise ResultsError("Every run must be an object")
        replicate = _positive_int(run.get("replicate"), "run.replicate")
        if replicate in seen_replicates or replicate > expected_n:
            raise ResultsError("Replicate ids must be unique and consecutive from 1")
        seen_replicates.add(replicate)
        if run.get("study_id", study_id) != study_id or run.get("reference", reference) != reference or run.get("comparison_type", "point_mutation") != "point_mutation":
            raise ResultsError("Run scope differs from the study/reference; cross-backbone pooling is forbidden")
        directory = _relative_path(root, run.get("directory"), "run.directory")
        diff_path = _relative_path(root, run.get("diff_file"), "run.diff_file")
        if not diff_path.is_relative_to(directory) or diff_path in seen_files:
            raise ResultsError("Each replicate must have its own difference output inside its directory")
        seen_files.add(diff_path)
        try:
            diff_digest = hashlib.sha256(diff_path.read_bytes()).hexdigest()
        except OSError as exc:
            raise ResultsError(f"Missing or unreadable replicate output {diff_path}: {exc}") from exc
        recorded_digest = run.get("diff_sha256")
        if not isinstance(recorded_digest, str) or len(recorded_digest) != 64 or any(char not in "0123456789abcdef" for char in recorded_digest):
            raise ResultsError(f"Replicate {replicate}: diff_sha256 must record the execution output digest")
        if recorded_digest != diff_digest:
            raise ResultsError(f"Replicate {replicate}: output SHA-256 does not match the recorded execution output")
        mapping = run.get("variants")
        if not isinstance(mapping, list) or len(mapping) != len(variants):
            raise ResultsError(f"Replicate {replicate} has an incomplete variant mapping")
        mapped: dict[int, dict[str, Any]] = {}
        for item in mapping:
            if not isinstance(item, dict):
                raise ResultsError("Run variant mappings must be objects")
            index = _positive_int(item.get("foldx_index"), "run.variant.foldx_index")
            if index not in variant_by_index or index in mapped or item.get("variant_id") != variant_by_index[index]["id"]:
                raise ResultsError(f"Replicate {replicate} has an inconsistent variant/index mapping")
            _text(item.get("expected_model_id"), "run.variant.expected_model_id")
            mapped[index] = item
        rows = parse_diff(diff_path)
        if len(rows) != len(variants):
            raise ResultsError(f"Replicate {replicate}: expected {len(variants)} difference rows, got {len(rows)}")
        if len({_model_key(row["model_id"]) for row in rows}) != len(rows):
            raise ResultsError(f"Replicate {replicate}: model ids collide after .pdb normalization")
        for index, row in enumerate(rows, 1):
            item = mapped[index]
            if _model_key(row["model_id"]) != _model_key(item["expected_model_id"]):
                raise ResultsError(f"Replicate {replicate}, FoldX position {index}: model id {row['model_id']!r} does not match {item['expected_model_id']!r}")
            measurements[item["variant_id"]].append({"replicate": replicate, "model_id": row["model_id"], "ddg_kcal_mol": row["ddg_kcal_mol"], "terms": row["terms"], "source_file": str(diff_path.relative_to(root))})
        provenance_files.append({"path": str(diff_path.relative_to(root)), "sha256": diff_digest, "replicate": replicate})
    _verify_step_artifacts(root, manifest.get("steps"), runs)
    summaries = []
    for index in sorted(variant_by_index):
        variant = variant_by_index[index]
        observed = sorted(measurements[variant["id"]], key=lambda item: item["replicate"])
        values = [item["ddg_kcal_mol"] for item in observed]
        summaries.append({"variant_id": variant["id"], "mutation_line": variant["mutation_line"], "foldx_index": index, "n": len(values), "mean_ddg_kcal_mol": statistics.mean(values), "sd_ddg_kcal_mol": statistics.stdev(values) if len(values) > 1 else None, "min_ddg_kcal_mol": min(values), "max_ddg_kcal_mol": max(values), "measurements": observed})
    return {"schema_version": 1, "study_id": study_id, "comparison_type": "point_mutation", "synthetic": manifest["synthetic"], "status": "completed", "reference": reference, "foldx": foldx, "variants": summaries, "uncertainty_label": "FoldX model sampling variability; not experimental uncertainty or a confidence interval", "sign_convention": "ΔΔG = G(mutant) − G(WT), kcal/mol; positive values predict destabilization", "provenance": {"manifest": "run_manifest.json", "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(), "source_files": provenance_files, "run_manifest": manifest}, "warnings": ["Replicates sample the FoldX modeling procedure; they do not establish accuracy or experimental reproducibility."]}
