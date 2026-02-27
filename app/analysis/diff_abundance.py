"""
MicrobiomeDash — Differential abundance analysis via R tools.
"""
import itertools
import json
import tempfile
import threading
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from biom import load_table as _load_biom

from app.analysis.r_runner import prepare_da_inputs, run_r_script
from app.analysis.shared import biom_to_count_df
from app.analysis.taxonomy import LEVEL_MAP, aggregate_counts_by_level

_TAX_RANKS = [k for k, v in sorted(LEVEL_MAP.items(), key=lambda x: x[1]) if v >= 0]


def _attach_taxonomy(df: pd.DataFrame, biom_path: str) -> pd.DataFrame:
    """Join full taxonomy lineage to DA results by feature (ASV) ID."""
    table = _load_biom(biom_path)
    tax_rows = []
    for obs_id in table.ids(axis="observation"):
        md = table.metadata(obs_id, axis="observation")
        if md and "taxonomy" in md:
            ranks = md["taxonomy"]
            if isinstance(ranks, str):
                ranks = [r.strip() for r in ranks.split(";")]
            row = {"feature": obs_id}
            for i, rank_name in enumerate(_TAX_RANKS):
                row[rank_name] = ranks[i] if i < len(ranks) and ranks[i] else ""
            tax_rows.append(row)

    if not tax_rows:
        return df

    tax_df = pd.DataFrame(tax_rows)
    return df.merge(tax_df, on="feature", how="left")


TOOL_SCRIPTS = {
    "ancombc": "run_ancombc.R",
    "aldex2": "run_aldex2.R",
    "deseq2": "run_deseq2.R",
    "linda": "run_linda.R",
    "maaslin2": "run_maaslin2.R",
}

TOOL_LABELS = {
    "ancombc": "ANCOM-BC2",
    "aldex2": "ALDEx2",
    "deseq2": "DESeq2",
    "linda": "LinDA",
    "maaslin2": "MaAsLin2",
}


def run_da_tool(
    tool: str,
    biom_path: str,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    matched_samples: list[str],
    group_col: str,
    ref_group: str,
    test_group: str,
    level: str = "ASV",
    threads: int | None = None,
    on_log=None,
) -> pd.DataFrame:
    """Run a differential abundance tool and return results DataFrame.

    Args:
        on_log: optional callback invoked with each stderr line from R for
                real-time log streaming.
    """
    script = TOOL_SCRIPTS.get(tool)
    if not script:
        raise ValueError(f"Unknown DA tool: {tool}")

    if level and level != "ASV":
        count_df = aggregate_counts_by_level(biom_path, level)
    else:
        count_df = biom_to_count_df(biom_path)
    output_dir = tempfile.mkdtemp()
    output_path = str(Path(output_dir) / "da_results.tsv")

    counts_path, meta_path = prepare_da_inputs(
        count_df, meta_df, sample_id_col, matched_samples,
        group_col, ref_group, test_group, output_dir,
    )

    script_args = {
        "counts": counts_path,
        "metadata": meta_path,
        "group_col": group_col,
        "ref_group": ref_group,
        "test_group": test_group,
        "output": output_path,
    }

    # Pass threads to tools that support parallelization
    # (ANCOM-BC2, ALDEx2, DESeq2, MaAsLin2 — not LinDA)
    if tool in ("ancombc", "aldex2", "deseq2", "maaslin2"):
        import os
        script_args["threads"] = threads or max(1, (os.cpu_count() or 1) // 2)

    result = run_r_script(script, script_args, on_line=on_log)

    if not result["success"]:
        raise RuntimeError(
            f"{TOOL_LABELS.get(tool, tool)} failed: {result.get('error', 'unknown')}"
        )

    if not Path(output_path).exists():
        raise RuntimeError(f"No output file produced by {TOOL_LABELS.get(tool, tool)}")

    df = pd.read_csv(output_path, sep="\t")

    # Attach full taxonomy lineage when running at ASV level
    if (not level or level == "ASV") and "feature" in df.columns:
        df = _attach_taxonomy(df, biom_path)

    return df


def build_volcano(
    results_df: pd.DataFrame,
    tool_name: str,
    q_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
) -> go.Figure:
    """Build a volcano plot from DA results."""
    df = results_df.copy()

    if "qvalue" not in df.columns or "log2fc" not in df.columns:
        return go.Figure().update_layout(title="No results to plot")

    df["neg_log10_q"] = -np.log10(df["qvalue"].clip(lower=1e-300))

    # Classify significance
    df["sig"] = "Not significant"
    sig_mask = (df["qvalue"] < q_thresh) & (df["log2fc"].abs() > lfc_thresh)
    df.loc[sig_mask & (df["log2fc"] > 0), "sig"] = "Up"
    df.loc[sig_mask & (df["log2fc"] < 0), "sig"] = "Down"

    colors = {"Not significant": "#636363", "Up": "#e74c3c", "Down": "#3498db"}

    fig = go.Figure()
    for sig_type in ["Not significant", "Up", "Down"]:
        sub = df[df["sig"] == sig_type]
        if sub.empty:
            continue
        fig.add_trace(go.Scatter(
            x=sub["log2fc"],
            y=sub["neg_log10_q"],
            mode="markers",
            name=sig_type,
            marker=dict(color=colors[sig_type], size=6, opacity=0.7),
            text=sub["feature"],
            hovertemplate="<b>%{text}</b><br>log2FC: %{x:.2f}<br>-log10(q): %{y:.2f}<extra></extra>",
        ))

    # Threshold lines
    fig.add_hline(y=-np.log10(q_thresh), line_dash="dash", line_color="gray", opacity=0.5)
    fig.add_vline(x=lfc_thresh, line_dash="dash", line_color="gray", opacity=0.5)
    fig.add_vline(x=-lfc_thresh, line_dash="dash", line_color="gray", opacity=0.5)

    fig.update_layout(
        title=f"Volcano Plot — {tool_name}",
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(q-value)",
        template="plotly_dark",
        height=500,
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
    )
    return fig


# ── All-pairwise background execution ───────────────────────────────────────


def _da_pairwise_progress_path(job_id: str) -> Path:
    return Path(f"/tmp/da_pairwise_{job_id}.json")


def _write_da_progress(job_id: str, data: dict) -> None:
    path = _da_pairwise_progress_path(job_id)
    path.write_text(json.dumps(data))


def read_da_pairwise_progress(job_id: str) -> dict | None:
    path = _da_pairwise_progress_path(job_id)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def cancel_da_pairwise(job_id: str) -> None:
    prog = read_da_pairwise_progress(job_id)
    if prog is None:
        prog = {}
    prog["cancel"] = True
    prog["status"] = "cancelled"
    _write_da_progress(job_id, prog)


def run_pairwise_da_background(
    tool: str,
    biom_path: str,
    meta_df: pd.DataFrame,
    sid_col: str,
    matched_samples: list[str],
    group_col: str,
    groups: list[str],
    job_id: str,
    level: str = "ASV",
    threads: int | None = None,
) -> None:
    """Spawn a background thread that runs DA for all C(n,2) group pairs."""

    def _worker():
        pairs = list(itertools.combinations(sorted(groups), 2))
        total = len(pairs)
        tool_label = TOOL_LABELS.get(tool, tool)
        progress = {
            "status": "running",
            "completed": 0,
            "total": total,
            "tool": tool,
            "tool_label": tool_label,
            "log": [],
            "results": [],
        }
        _write_da_progress(job_id, progress)

        for i, (ref, test) in enumerate(pairs):
            # Check cancellation
            current = read_da_pairwise_progress(job_id)
            if current and current.get("cancel"):
                progress["status"] = "cancelled"
                _write_da_progress(job_id, progress)
                return

            label = f"{ref} vs {test}"
            progress["log"].append(
                f"[{i+1}/{total}] Running {tool_label}: {label}..."
            )
            _write_da_progress(job_id, progress)

            def _on_log(line, _prog=progress, _job=job_id):
                _prog["log"].append(f"  {line}")
                _write_da_progress(_job, _prog)

            try:
                df = run_da_tool(
                    tool, biom_path, meta_df, sid_col, matched_samples,
                    group_col, ref, test, level=level, threads=threads,
                    on_log=_on_log,
                )
                df["comparison"] = label
                df["ref_group"] = ref
                df["test_group"] = test

                n_feat = len(df)
                n_sig = int((df["qvalue"] < 0.05).sum()) if "qvalue" in df.columns else 0
                progress["log"].append(
                    f"[{i+1}/{total}] {label} — {n_feat} features, {n_sig} significant"
                )
                progress["results"].extend(df.to_dict(orient="records"))
            except Exception as e:
                progress["log"].append(f"[{i+1}/{total}] {label} — ERROR: {e}")

            progress["completed"] = i + 1
            _write_da_progress(job_id, progress)

        progress["status"] = "complete"
        _write_da_progress(job_id, progress)

    thread = threading.Thread(target=_worker, daemon=True)
    thread.start()
