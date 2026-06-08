#!/usr/bin/env python3
"""
Node 02 — Visualize Report

Generates an interactive HTML solubility report containing:
  - Summary table (best solvent per molecule, top-N heatmap)
  - Colour-coded heatmap: molecules × top-N solvents by mean logS
  - Per-molecule horizontal bar charts ranked by logS
"""

import argparse
import json
import os
import sys

import pandas as pd
import numpy as np


# ── Solubility category thresholds (log mol/L) ──────────────────────────────
def _logS_label(logS):
    if logS is None or np.isnan(logS):
        return ("—", "#cccccc")
    if logS >= 0:
        return ("≥1 mol/L", "#1a7f1a")
    if logS >= -2:
        return ("10–1000 mmol/L", "#4caf50")
    if logS >= -4:
        return ("0.1–10 mmol/L", "#ff9800")
    if logS >= -6:
        return ("1–100 µmol/L", "#f44336")
    return ("<1 µmol/L", "#9c27b0")


def _color_scale(val, vmin, vmax):
    """Map a logS value to a green–yellow–red hex colour."""
    if val is None or np.isnan(val):
        return "#eeeeee"
    frac = (val - vmin) / (vmax - vmin + 1e-9)
    frac = max(0.0, min(1.0, frac))
    r = int(255 * (1 - frac))
    g = int(200 * frac + 50)
    b = 40
    return f"rgb({r},{g},{b})"


def _bar_chart_html(mol_name, mol_df):
    """Return an inline Plotly bar-chart HTML block for one molecule."""
    mol_df = mol_df.dropna(subset=["logS"]).sort_values("logS", ascending=False)
    if mol_df.empty:
        return "<p>No predictions available.</p>"

    solvents = mol_df["solvent"].tolist()
    logS     = mol_df["logS"].tolist()
    stdev    = mol_df["logS_stdev"].fillna(0).tolist()
    temp_c   = mol_df["temperature_c"].iloc[0]

    vmin, vmax = min(logS), max(logS)
    colors = [_color_scale(v, vmin, vmax) for v in logS]

    bars = []
    for sv, ls, sd, col in zip(solvents, logS, stdev, colors):
        bars.append(
            f'{{"x":{ls:.4f},"y":"{sv}","error_x":{sd:.4f},"color":"{col}"}}'
        )
    bars_json = "[" + ",".join(bars) + "]"

    chart_id = mol_name.replace(" ", "_").replace("/", "_")
    return f"""
<div class="chart-block">
  <h3>{mol_name} <span class="temp-badge">{temp_c:.0f} °C</span></h3>
  <div id="chart_{chart_id}" style="width:100%;height:{max(300, len(solvents)*22)}px;"></div>
  <script>
  (function(){{
    var data = {bars_json};
    var trace = {{
      type: 'bar', orientation: 'h',
      x: data.map(d=>d.x),
      y: data.map(d=>d.y),
      error_x: {{type:'data', array:data.map(d=>d.error_x), visible:true,
                 color:'rgba(80,80,80,0.5)', thickness:1.5, width:4}},
      marker: {{color: data.map(d=>d.color)}},
      hovertemplate: '<b>%{{y}}</b><br>logS = %{{x:.4f}}<extra></extra>',
    }};
    var layout = {{
      margin: {{l:160, r:40, t:20, b:50}},
      xaxis: {{title:'log(S) [mol/L]', zeroline:true, zerolinewidth:2,
               zerolinecolor:'rgba(0,180,0,0.5)'}},
      yaxis: {{automargin:true}},
      shapes: [
        {{type:'line', x0:0,   x1:0,   y0:-0.5, y1:data.length-0.5,
          line:{{color:'rgba(0,180,0,0.35)',  width:2, dash:'dot'}}}},
        {{type:'line', x0:-3,  x1:-3,  y0:-0.5, y1:data.length-0.5,
          line:{{color:'rgba(255,140,0,0.35)', width:2, dash:'dot'}}}},
        {{type:'line', x0:-6,  x1:-6,  y0:-0.5, y1:data.length-0.5,
          line:{{color:'rgba(200,0,0,0.35)',   width:2, dash:'dot'}}}},
      ],
      paper_bgcolor:'rgba(0,0,0,0)',
      plot_bgcolor:'rgba(248,248,255,1)',
    }};
    Plotly.newPlot('chart_{chart_id}', [trace], layout, {{responsive:true, displayModeBar:false}});
  }})();
  </script>
</div>
"""


def _heatmap_html(pivot, vmin, vmax):
    """Render an HTML table heatmap for molecule × solvent."""
    solvents = list(pivot.columns)
    molecules = list(pivot.index)

    header_cells = "".join(
        f'<th title="{sv}">{sv[:18]}{"…" if len(sv)>18 else ""}</th>'
        for sv in solvents
    )
    rows_html = ""
    for mol in molecules:
        cells = ""
        for sv in solvents:
            val = pivot.loc[mol, sv]
            bg = _color_scale(val, vmin, vmax)
            txt = f"{val:.2f}" if not np.isnan(val) else "—"
            cells += f'<td style="background:{bg};color:#fff;font-size:11px;" title="{sv}: logS={txt}">{txt}</td>'
        rows_html += f"<tr><td class='mol-name'>{mol}</td>{cells}</tr>"

    return f"""
<div class="heatmap-wrap">
  <table class="heatmap">
    <thead><tr><th>Molecule</th>{header_cells}</tr></thead>
    <tbody>{rows_html}</tbody>
  </table>
</div>
"""


def _summary_table_html(df):
    """Best solvent per molecule summary."""
    best = (
        df.dropna(subset=["logS"])
          .sort_values("logS", ascending=False)
          .groupby("name", sort=False)
          .first()
          .reset_index()
    )
    rows = ""
    for _, row in best.iterrows():
        lbl, _ = _logS_label(row["logS"])
        rows += (
            f"<tr>"
            f"<td>{row['name']}</td>"
            f"<td><code>{row['smiles']}</code></td>"
            f"<td>{row['solvent']}</td>"
            f"<td style='font-family:monospace'>{row['logS']:+.3f}</td>"
            f"<td style='font-family:monospace'>±{row['logS_stdev']:.3f}</td>"
            f"<td>{lbl}</td>"
            f"</tr>"
        )
    return f"""
<table class='summary-table'>
  <thead>
    <tr>
      <th>Molecule</th><th>SMILES</th><th>Best Solvent</th>
      <th>logS</th><th>Stdev</th><th>Category</th>
    </tr>
  </thead>
  <tbody>{rows}</tbody>
</table>
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top_n_solvents", type=int, default=20)
    args = parser.parse_args()

    # --- Load input ---
    input_path = "./inputs/solubility_predictions.csv"
    if not os.path.exists(input_path):
        print("ERROR: solubility_predictions.csv not found in inputs/", flush=True)
        sys.exit(1)

    df = pd.read_csv(input_path)
    molecules = df["name"].unique().tolist()
    temp_c = df["temperature_c"].iloc[0] if "temperature_c" in df.columns else 25.0
    print(f"Building report: {len(molecules)} molecules, {df['solvent'].nunique()} solvents.", flush=True)

    # --- Select top-N solvents by mean logS across all molecules ---
    mean_logS = df.groupby("solvent")["logS"].mean().sort_values(ascending=False)
    top_solvents = mean_logS.head(args.top_n_solvents).index.tolist()

    # --- Heatmap pivot ---
    pivot = df[df["solvent"].isin(top_solvents)].pivot_table(
        index="name", columns="solvent", values="logS", aggfunc="mean"
    ).reindex(columns=top_solvents)
    vmin = float(df["logS"].min())
    vmax = float(df["logS"].max())

    # --- Per-molecule bar charts ---
    charts_html = ""
    for mol in molecules:
        mol_df = df[df["name"] == mol].copy()
        charts_html += _bar_chart_html(mol, mol_df)

    # --- Assemble HTML ---
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Solubility Report — fastsolv</title>
  <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            margin: 0; padding: 24px 32px; background: #f7f8fc; color: #222; }}
    h1   {{ font-size: 1.6rem; color: #1a3a6b; border-bottom: 2px solid #1a3a6b;
            padding-bottom: 8px; }}
    h2   {{ font-size: 1.2rem; color: #1a3a6b; margin-top: 32px; }}
    h3   {{ font-size: 1rem; color: #333; margin: 20px 0 4px; }}
    .meta {{ font-size: 13px; color: #666; margin-bottom: 20px; }}
    .temp-badge {{ background: #1a3a6b; color: #fff; border-radius: 4px;
                   padding: 2px 6px; font-size: 11px; font-weight: normal; }}
    /* Summary table */
    .summary-table {{ border-collapse: collapse; width: 100%; font-size: 13px; margin-bottom: 24px; }}
    .summary-table th {{ background: #1a3a6b; color: #fff; padding: 6px 10px; text-align: left; }}
    .summary-table td {{ padding: 5px 10px; border-bottom: 1px solid #e0e0e0; }}
    .summary-table tr:nth-child(even) {{ background: #f0f4ff; }}
    code {{ font-size: 11px; color: #555; word-break: break-all; }}
    /* Heatmap */
    .heatmap-wrap {{ overflow-x: auto; margin-bottom: 32px; }}
    .heatmap {{ border-collapse: collapse; font-size: 11px; white-space: nowrap; }}
    .heatmap th {{ background: #1a3a6b; color: #fff; padding: 4px 6px;
                   writing-mode: vertical-rl; text-orientation: mixed;
                   transform: rotate(180deg); max-height: 120px; }}
    .heatmap th:first-child {{ writing-mode: horizontal-tb; transform: none;
                                min-width: 120px; }}
    .heatmap td {{ padding: 4px 8px; border: 1px solid rgba(255,255,255,0.2); text-align: center; }}
    .mol-name {{ background: #f0f4ff !important; color: #222 !important;
                 font-weight: 600; text-align: left !important; padding: 4px 10px !important; }}
    /* Legend */
    .legend {{ display: flex; gap: 16px; flex-wrap: wrap; margin: 8px 0 24px; font-size: 12px; }}
    .legend-item {{ display: flex; align-items: center; gap: 6px; }}
    .legend-swatch {{ width: 16px; height: 16px; border-radius: 3px; }}
    /* Charts */
    .chart-block {{ background: #fff; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.08);
                    padding: 16px 20px; margin-bottom: 24px; }}
  </style>
</head>
<body>
  <h1>🧪 Solubility Report</h1>
  <div class="meta">
    Generated by <strong>fastsolv</strong> (MIT, 2025) &nbsp;|&nbsp;
    <strong>{len(molecules)}</strong> molecules &nbsp;|&nbsp;
    <strong>{df['solvent'].nunique()}</strong> solvents &nbsp;|&nbsp;
    Temperature: <strong>{temp_c:.0f} °C</strong>
  </div>

  <div class="legend">
    <span style="font-weight:600">log(S) key:</span>
    <div class="legend-item"><div class="legend-swatch" style="background:#1a7f1a"></div>≥1 mol/L</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#4caf50"></div>10–1000 mmol/L</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#ff9800"></div>0.1–10 mmol/L</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#f44336"></div>1–100 µmol/L</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#9c27b0"></div>&lt;1 µmol/L</div>
  </div>

  <h2>Best Solvent per Molecule</h2>
  {_summary_table_html(df)}

  <h2>Heatmap — Top {args.top_n_solvents} Solvents by Mean Solubility</h2>
  {_heatmap_html(pivot, vmin, vmax)}

  <h2>Solubility Profiles</h2>
  {charts_html}

</body>
</html>
"""

    os.makedirs("./outputs", exist_ok=True)
    with open("./outputs/solubility_report.html", "w", encoding="utf-8") as f:
        f.write(html)
    print("Saved solubility_report.html.", flush=True)


if __name__ == "__main__":
    main()
