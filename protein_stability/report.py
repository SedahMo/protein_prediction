"""Dependency-free, self-contained stability reports and a synthetic demo."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import math
from pathlib import Path
from typing import Any

from .results import collect_results


COMPARISON_MATRIX = [
    ("point_mutation", "Point substitutions", "Same sequence length, domain order and reference backbone; explicit mutations", "BuildModel mutant − WT ΔΔG within this study; compare paired variants only"),
    ("linker", "Different linker", "Separate validated structural ensembles; inspect junction geometry and domain contacts", "Assess ensemble quality, developability and experiments; do not rank absolute FoldX energies across constructs"),
    ("domain_order", "VH/VL order change", "Model both orientations and linker contexts; evaluate interface confidence", "Use orientation-specific evidence and matched experiments; no direct BuildModel ΔΔG for rearrangement"),
    ("cross_protein", "Different proteins/backbones", "Independent sequence/structure provenance, residue coverage and structural confidence", "No universal stability ranking from FoldX total energy; use appropriate matched endpoints"),
]


def blocked_report(study_id: str, reason: str, comparison_type: str = "point_mutation") -> dict[str, Any]:
    """Return a reportable template when a licensed tool or structural input is absent."""
    return {"schema_version": 1, "study_id": study_id, "comparison_type": comparison_type, "synthetic": False, "status": "blocked", "variants": [], "warnings": [reason], "provenance": {}, "next_steps": ["Supply a licensed FoldX executable and validated, sequence-matched reference PDB before producing energies."]}


def _e(value: Any) -> str:
    return html.escape(str(value), quote=True)


def _num(value: Any) -> str:
    return "—" if value is None else f"{value:.3f}"


def _csv_safe(value: Any) -> Any:
    return "'" + value if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")) else value


def _chart(variants: list[dict[str, Any]]) -> str:
    if not variants:
        return ""
    width, left, right = 900, 220, 60
    height = 95 + 48 * len(variants)
    limit = max(abs(v["mean_ddg_kcal_mol"]) + (v.get("sd_ddg_kcal_mol") or 0) for v in variants)
    limit = max(0.5, limit * 1.2)
    scale = (width - left - right) / (2 * limit)
    zero = left + (width - left - right) / 2
    marks = []
    for fraction in (-1, -0.5, 0, 0.5, 1):
        val = limit * fraction
        x = zero + val * scale
        marks.append(f'<line x1="{x:.2f}" y1="28" x2="{x:.2f}" y2="{height-45}" stroke="{"#768176" if fraction == 0 else "#e1ded4"}"/><text x="{x:.2f}" y="{height-25}" text-anchor="middle">{val:.2f}</text>')
    for index, variant in enumerate(variants):
        y = 52 + index * 48
        mean = variant["mean_ddg_kcal_mol"]
        sd = variant.get("sd_ddg_kcal_mol") or 0
        x, x1, x2 = zero + mean * scale, zero + (mean - sd) * scale, zero + (mean + sd) * scale
        color = "#376b57" if mean < 0 else "#a35a35"
        marks.append(f'<text x="{left-18}" y="{y+5}" text-anchor="end">{_e(variant["variant_id"])}</text><line x1="{x1:.2f}" y1="{y}" x2="{x2:.2f}" y2="{y}" stroke="{color}" stroke-width="3"/><line x1="{x1:.2f}" y1="{y-5}" x2="{x1:.2f}" y2="{y+5}" stroke="{color}"/><line x1="{x2:.2f}" y1="{y-5}" x2="{x2:.2f}" y2="{y+5}" stroke="{color}"/><circle cx="{x:.2f}" cy="{y}" r="6" fill="{color}"/>')
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="FoldX point mutation mean difference with sample standard deviation"><title>Point mutation ΔΔG: mean ± sample SD</title><g font-family="system-ui,sans-serif" font-size="13" fill="#303b35">{"".join(marks)}<text x="{zero}" y="{height-3}" text-anchor="middle">ΔΔG (kcal/mol): positive predicts destabilization</text></g></svg>'


def write_report(data: dict[str, Any], outdir: str | Path) -> dict[str, str]:
    """Write JSON, CSV and an offline HTML report without creating a ranking.

    Non-point comparisons and blocked studies get an explicit input/status
    template, and never receive the point-mutation chart or ΔΔG summary table.
    """
    # Serialization first rejects nonfinite values before any artifacts are made.
    serialized = json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    outdir = Path(outdir)
    variants = data.get("variants", [])
    if not isinstance(variants, list):
        raise ValueError("Report variants must be a list")
    eligible = data.get("comparison_type") == "point_mutation" and data.get("status") in ("complete", "completed")
    if variants and not eligible:
        raise ValueError("Energy summaries are not eligible for blocked or non-point comparisons")
    for variant in variants:
        for key in ("mean_ddg_kcal_mol", "min_ddg_kcal_mol", "max_ddg_kcal_mol"):
            value = variant.get(key)
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value):
                raise ValueError(f"Invalid numeric summary: {key}")
        sd = variant.get("sd_ddg_kcal_mol")
        if sd is not None and (not isinstance(sd, (int, float)) or isinstance(sd, bool) or not math.isfinite(sd) or sd < 0):
            raise ValueError("Invalid sample standard deviation")
    outdir.mkdir(parents=True, exist_ok=True)
    json_path, csv_path, html_path = outdir / "results.json", outdir / "summary.csv", outdir / "report.html"
    json_path.write_text(serialized, encoding="utf-8")
    fields = ["study_id", "comparison_type", "synthetic", "variant_id", "mutation_line", "n", "mean_ddg_kcal_mol", "sd_ddg_kcal_mol", "min_ddg_kcal_mol", "max_ddg_kcal_mol"]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for variant in variants:
            row = {field: variant.get(field, "") for field in fields}
            row.update({"study_id": data.get("study_id", ""), "comparison_type": data.get("comparison_type", ""), "synthetic": data.get("synthetic", False)})
            writer.writerow({key: _csv_safe(value) for key, value in row.items()})
    synthetic = data.get("synthetic") is True
    if synthetic:
        banner = '<div class="banner synthetic"><strong>SYNTHETIC DEMONSTRATION · 人工演示数据</strong><br>Hand-written fixture values test this workflow. No biological stability calculation was performed.</div>'
    elif not eligible or not variants:
        banner = '<div class="banner"><strong>ANALYSIS TEMPLATE · 待输入 / 待计算</strong><br>No completed energy calculations are reported. 此页面用于记录输入、计算状态和下一步工作。</div>'
    else:
        banner = '<div class="banner"><strong>COMPUTATIONAL OUTPUT · 计算结果</strong><br>Provenance and run status below define which calculations were completed.</div>'
    warning_html = "".join(f"<li>{_e(warning)}</li>" for warning in data.get("warnings", []))
    next_html = "".join(f"<li>{_e(step)}</li>" for step in data.get("next_steps", []))
    matrix = "".join(f'<tr class="{"active" if data.get("comparison_type") == mode else ""}"><td>{_e(name)}<br><code>{mode}</code></td><td>{_e(inputs)}</td><td>{_e(interpretation)}</td></tr>' for mode, name, inputs, interpretation in COMPARISON_MATRIX)
    summary_rows = "".join(f'<tr><td>{_e(v["variant_id"])}</td><td><code>{_e(v.get("mutation_line", ""))}</code></td><td>{_e(v.get("n", ""))}</td><td>{_num(v["mean_ddg_kcal_mol"])}</td><td>{_num(v.get("sd_ddg_kcal_mol"))}</td><td>{_num(v["min_ddg_kcal_mol"])}</td><td>{_num(v["max_ddg_kcal_mol"])}</td></tr>' for v in variants)
    results_html = ('<section><h2>Within-reference point mutation comparison</h2><p>ΔΔG = G(mutant) − G(WT), kcal/mol. Positive values predict destabilization. Variants retain their manifest order.</p>' + _chart(variants) + '<div class="scroll"><table><thead><tr><th>Variant</th><th>Mutation input</th><th>n</th><th>Mean</th><th>Sample SD</th><th>Minimum</th><th>Maximum</th></tr></thead><tbody>' + summary_rows + '</tbody></table></div><p class="note">Mean ± sample SD describes FoldX model sampling variability. It is not experimental uncertainty or a confidence interval. SD is unavailable for a single repeat. This report does not convert ΔΔG into melting temperature.</p></section>' if eligible and variants else '<section><h2>Analysis status</h2><p>No eligible point mutation energies are reported. Supply the required inputs and complete the defined run before interpreting stability.</p></section>')
    measurement_html = []
    for variant in variants:
        measurement_rows = "".join(f'<tr><td>{_e(m.get("replicate", ""))}</td><td>{_e(m.get("model_id", ""))}</td><td>{_num(m["ddg_kcal_mol"])}</td><td><code>{_e(m.get("source_file", ""))}</code></td></tr>' for m in variant.get("measurements", []))
        measurement_html.append(f'<details><summary>{_e(variant["variant_id"])} — replicate trace</summary><div class="scroll"><table><thead><tr><th>Repeat</th><th>Model</th><th>ΔΔG (kcal/mol)</th><th>Source</th></tr></thead><tbody>{measurement_rows}</tbody></table></div></details>')
    provenance = {"reference": data.get("reference", {}), "foldx": data.get("foldx", {}), "provenance": data.get("provenance", {})}
    document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; img-src data:; base-uri 'none'; form-action 'none'"><title>{_e(data.get('study_id', 'Stability study'))} · Stability analysis</title>
<style>:root{{color-scheme:light}}*{{box-sizing:border-box}}body{{margin:0;background:#f4f0e6;color:#28362e;font:16px/1.65 system-ui,-apple-system,sans-serif}}main{{max-width:1180px;margin:40px auto;padding:0 28px 60px}}header{{border-top:5px solid #376b57;padding:26px 0 18px}}h1{{font-family:Georgia,serif;font-size:38px;line-height:1.18;margin:8px 0 16px}}h2{{font-family:Georgia,serif;font-size:25px;margin:0 0 15px}}.eyebrow{{font-size:12px;letter-spacing:.16em;text-transform:uppercase;color:#637468}}.meta{{display:flex;gap:12px;flex-wrap:wrap}}.chip{{padding:4px 10px;border:1px solid #c6cec2;border-radius:4px;font-size:13px}}section{{background:#fffcf4;border:1px solid #dedccd;padding:24px;margin:22px 0}}.banner{{padding:17px 20px;border:1px solid #a4b6a1;background:#e6eddf;margin:14px 0}}.synthetic{{background:#fff0ce;border:2px solid #ac6a16;color:#73430e}}table{{border-collapse:collapse;width:100%;font-size:14px;text-align:left}}td,th{{padding:12px 10px;border-bottom:1px solid #dedccd;vertical-align:top}}th{{background:#eeeee3;font-size:12px;letter-spacing:.03em}}tr.active td{{background:#edf2e7}}.scroll{{overflow:auto}}.note{{font-size:13px;color:#687267}}code,pre{{font-family:ui-monospace,SFMono-Regular,monospace;font-size:12px}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#efeee4;padding:18px}}svg{{display:block;width:100%;min-width:540px}}details{{padding:12px 0;border-top:1px solid #dedccd}}summary{{cursor:pointer}}footer{{font-size:12px;color:#687267}}@media print{{body{{background:white}}main{{margin:0;padding:0}}section{{break-inside:avoid}}details>div{{display:block}}}}@media(max-width:600px){{main{{padding:0 15px}}h1{{font-size:30px}}section{{padding:16px}}}}</style></head>
<body><main><header><div class="eyebrow">Antibody engineering · reproducible computation</div><h1>{_e(data.get('study_id', 'Stability study'))}</h1><div class="meta"><span class="chip">Status: {_e(data.get('status', 'unspecified'))}</span><span class="chip">Comparison: {_e(data.get('comparison_type', 'unspecified'))}</span></div></header>{banner}{results_html}
<section><h2>Comparison rules</h2><div class="scroll"><table><thead><tr><th>Comparison</th><th>Required basis</th><th>Permitted conclusion</th></tr></thead><tbody>{matrix}</tbody></table></div></section>
<section><h2>Interpretation and remaining work</h2><ul>{warning_html or '<li>Review sequence–structure matching, residue coverage and model quality with the recorded inputs.</li>'}</ul>{'<ol>'+next_html+'</ol>' if next_html else ''}</section>
{'<section><h2>Replicate outputs</h2>'+''.join(measurement_html)+'</section>' if variants else ''}
<section><h2>Traceability</h2><p>Energy data, manifest hashes, reference hashes and settings are exported in <code>results.json</code>; numeric summaries are exported in <code>summary.csv</code>.</p><details><summary>Reference, settings and source manifest</summary><pre>{_e(json.dumps(provenance, ensure_ascii=False, indent=2, allow_nan=False))}</pre></details></section><footer>Offline report · no external assets or network requests · generated by protein_stability</footer></main></body></html>'''
    html_path.write_text(document, encoding="utf-8")
    return {"results_json": str(json_path), "summary_csv": str(csv_path), "report_html": str(html_path)}


def make_demo(outdir: str | Path) -> dict[str, str]:
    """Generate conspicuously synthetic fixtures, collect them and write a report."""
    outdir = Path(outdir)
    fixture_dir = outdir / "synthetic_run"
    fixture_dir.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(b"SYNTHETIC FIXTURE: no structure or biological energy").hexdigest()
    variants = [{"id": "synthetic_variant_A", "foldx_index": 1, "mutation_line": "AA1V;"}, {"id": "synthetic_variant_B", "foldx_index": 2, "mutation_line": "GA2S;"}]
    values = [[-0.4, -0.2, -0.3, -0.1, -0.5], [0.2, 0.4, 0.3, 0.5, 0.6]]
    runs = []
    for index in range(5):
        relative_dir = f"replicates/run_{index+1:03d}"
        directory = fixture_dir / relative_dir
        directory.mkdir(parents=True, exist_ok=True)
        diff_file = f"{relative_dir}/Dif_demo_Repair.fxout"
        (fixture_dir / diff_file).write_text("SYNTHETIC TEST FIXTURE — NO BIOLOGICAL CALCULATION\nPdb\ttotal energy\tBackbone Hbond\tSidechain Hbond\n" + "".join(f"demo_Repair_{variant_index+1}.pdb\t{values[variant_index][index]}\t0.0\t0.0\n" for variant_index in range(2)), encoding="utf-8")
        runs.append({"replicate": index + 1, "directory": relative_dir, "diff_file": diff_file, "diff_sha256": hashlib.sha256((fixture_dir / diff_file).read_bytes()).hexdigest(), "variants": [{"variant_id": variant["id"], "foldx_index": variant["foldx_index"], "expected_model_id": f"demo_Repair_{variant['foldx_index']}.pdb"} for variant in variants]})
    manifest = {"schema_version": 1, "study_id": "SYNTHETIC DEMO — pipeline verification only", "comparison_type": "point_mutation", "synthetic": True, "status": "completed", "reference": {"pdb_sha256": digest, "fasta_sha256": digest, "prepared_pdb_sha256": digest}, "foldx": {"runs": 5, "temperature": 298.0, "pH": 7.0, "ion_strength": 0.05, "version": "NOT RUN: hand-written synthetic fixture", "seed_policy": "Not applicable: synthetic values"}, "variants": variants, "runs": runs, "notes": ["Hashes identify a synthetic fixture label, not an existing PDB or FASTA.", "Mutation labels and all numerical energies are artificial test data."]}
    manifest_path = fixture_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    paths = write_report(collect_results(fixture_dir), outdir)
    paths["run_manifest"] = str(manifest_path)
    return paths
