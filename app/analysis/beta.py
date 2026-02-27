"""
MicrobiomeDash — Beta diversity computation: distance matrices, ordination, PERMANOVA.
"""
import json
import tempfile
import threading
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from skbio import DistanceMatrix
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from statsmodels.stats.multitest import multipletests

from app.analysis.r_runner import run_r_script


def compute_distance(count_df: pd.DataFrame, metric: str) -> DistanceMatrix:
    """Compute a beta diversity distance matrix.

    Args:
        count_df: Features (rows) x samples (columns) integer count DataFrame.
        metric: 'braycurtis' or 'jaccard'.
    """
    counts = count_df.T  # samples as rows
    sample_ids = counts.index.tolist()
    dm = beta_diversity(metric, counts.values.astype(int), ids=sample_ids)
    return dm


def run_pcoa(dm: DistanceMatrix) -> tuple[pd.DataFrame, dict]:
    """Run PCoA ordination.

    Returns (coordinates DataFrame with PC1/PC2 columns, proportion_explained dict).
    """
    result = pcoa(dm)
    coords = result.samples[["PC1", "PC2"]].copy()
    coords.index.name = "sample_id"
    prop_explained = {
        f"PC{i+1}": float(v) for i, v in enumerate(result.proportion_explained[:2])
    }
    return coords, prop_explained


def compute_confidence_ellipse(
    coords: pd.DataFrame,
    group_values: pd.Series,
    confidence: float = 0.95,
) -> list[dict]:
    """Compute 95% confidence ellipses for each group on a 2D ordination.

    Args:
        coords: DataFrame with 'Axis1', 'Axis2' columns indexed by sample ID.
        group_values: Series mapping sample ID → group label (same index as coords).
        confidence: Confidence level (default 0.95).

    Returns:
        List of {"group": str, "x": list[float], "y": list[float]} for each group
        with ≥3 samples.
    """
    scale = chi2.ppf(confidence, df=2)
    theta = np.linspace(0, 2 * np.pi, 100)
    circle = np.column_stack([np.cos(theta), np.sin(theta)])

    ellipses = []
    for grp in sorted(group_values.unique()):
        mask = group_values == grp
        pts = coords.loc[mask, ["Axis1", "Axis2"]].values
        if len(pts) < 3:
            continue
        center = pts.mean(axis=0)
        cov = np.cov(pts, rowvar=False)
        eigvals, eigvecs = np.linalg.eigh(cov)
        # Scale eigenvectors by sqrt(eigenvalue * chi2_scale)
        radii = np.sqrt(eigvals * scale)
        ellipse_pts = circle * radii  # (100, 2)
        ellipse_pts = ellipse_pts @ eigvecs.T + center
        ellipses.append({
            "group": grp,
            "x": ellipse_pts[:, 0].tolist(),
            "y": ellipse_pts[:, 1].tolist(),
        })
    return ellipses


def run_nmds(dm: DistanceMatrix, output_dir: str | None = None) -> tuple[pd.DataFrame, float]:
    """Run NMDS ordination via R's vegan::metaMDS.

    Returns (coordinates DataFrame with NMDS1/NMDS2, stress value).
    """
    if output_dir is None:
        output_dir = tempfile.mkdtemp()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write distance matrix as TSV
    dm_path = str(output_dir / "distance_matrix.tsv")
    dm_df = dm.to_data_frame()
    dm_df.to_csv(dm_path, sep="\t")

    output_path = str(output_dir / "nmds_coords.tsv")

    result = run_r_script("run_nmds.R", {
        "dist_matrix": dm_path,
        "output": output_path,
        "k": 2,
    })

    if not result["success"]:
        raise RuntimeError(f"NMDS failed: {result.get('error', 'unknown')}")

    coords = pd.read_csv(output_path, sep="\t", index_col=0)
    stress = result.get("status", {}).get("stress", 0.0)
    return coords, float(stress)


def run_permanova_global(
    dm: DistanceMatrix,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
    n_permutations: int = 999,
) -> dict | None:
    """Run a single global PERMANOVA test.

    Returns {n_groups, n_samples, pseudo_F, pvalue} or None on failure.
    """
    dm_ids = dm.ids
    meta = meta_df.set_index(meta_df[sample_id_col].astype(str))
    meta = meta.loc[meta.index.intersection(dm_ids)]

    groups = meta[group_col].unique().tolist()
    if len(groups) < 2:
        return None

    common = [s for s in dm_ids if s in meta.index]
    if len(common) < 3:
        return None

    sub_dm = dm.filter(common)
    grouping = meta.loc[common, group_col]

    try:
        res = permanova(sub_dm, grouping, permutations=n_permutations)
        return {
            "n_groups": len(groups),
            "n_samples": len(common),
            "pseudo_F": float(res["test statistic"]),
            "pvalue": float(res["p-value"]),
        }
    except Exception:
        return None


def run_permanova_pairwise(
    dm: DistanceMatrix,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
    n_permutations: int = 999,
) -> pd.DataFrame:
    """Run pairwise PERMANOVA between all group pairs with BH correction.

    Returns DataFrame with group1, group2, n, pseudo_F, pvalue, qvalue.
    """
    # Align metadata to DM sample order
    dm_ids = dm.ids
    meta = meta_df.set_index(meta_df[sample_id_col].astype(str))
    meta = meta.loc[meta.index.intersection(dm_ids)]

    groups = meta[group_col].unique().tolist()
    if len(groups) < 2:
        return pd.DataFrame(columns=["group1", "group2", "n", "pseudo_F", "pvalue", "qvalue"])

    results = []
    p_values = []

    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            g1, g2 = groups[i], groups[j]
            subset_ids = meta[meta[group_col].isin([g1, g2])].index.tolist()
            subset_ids = [s for s in subset_ids if s in dm_ids]
            if len(subset_ids) < 3:
                continue

            sub_dm = dm.filter(subset_ids)
            grouping = meta.loc[subset_ids, group_col]

            try:
                res = permanova(sub_dm, grouping, permutations=n_permutations)
                results.append({
                    "group1": g1,
                    "group2": g2,
                    "n": len(subset_ids),
                    "pseudo_F": float(res["test statistic"]),
                    "pvalue": float(res["p-value"]),
                })
                p_values.append(float(res["p-value"]))
            except Exception:
                continue

    if results and p_values:
        _, q_values, _, _ = multipletests(p_values, method="fdr_bh")
        for row, q in zip(results, q_values):
            row["qvalue"] = float(q)

    return pd.DataFrame(results)


def _permanova_progress_path(job_id: str) -> Path:
    return Path(f"/tmp/bd_permanova_{job_id}.json")


def _write_progress(job_id: str, data: dict):
    path = _permanova_progress_path(job_id)
    path.write_text(json.dumps(data))


def run_permanova_pairwise_background(
    dm: DistanceMatrix,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
    n_permutations: int,
    job_id: str,
):
    """Run pairwise PERMANOVA in a background thread with progress tracking.

    Writes progress JSON to /tmp/bd_permanova_{job_id}.json after each pair.
    """
    def _worker():
        try:
            dm_ids = dm.ids
            meta = meta_df.set_index(meta_df[sample_id_col].astype(str))
            meta = meta.loc[meta.index.intersection(dm_ids)]
            groups = sorted(meta[group_col].unique().tolist())

            pairs = [(groups[i], groups[j])
                     for i in range(len(groups))
                     for j in range(i + 1, len(groups))]
            total = len(pairs)

            results = []
            p_values = []
            log = []

            _write_progress(job_id, {
                "status": "running", "completed": 0, "total": total,
                "log": log, "results": [],
            })

            for idx, (g1, g2) in enumerate(pairs):
                subset_ids = meta[meta[group_col].isin([g1, g2])].index.tolist()
                subset_ids = [s for s in subset_ids if s in dm_ids]

                if len(subset_ids) < 3:
                    log.append(f"[{idx+1}/{total}] {g1} vs {g2} — skipped (<3 samples)")
                    _write_progress(job_id, {
                        "status": "running", "completed": idx + 1, "total": total,
                        "log": log, "results": results,
                    })
                    continue

                try:
                    sub_dm = dm.filter(subset_ids)
                    grouping = meta.loc[subset_ids, group_col]
                    res = permanova(sub_dm, grouping, permutations=n_permutations)
                    row = {
                        "group1": g1, "group2": g2,
                        "n": len(subset_ids),
                        "pseudo_F": float(res["test statistic"]),
                        "pvalue": float(res["p-value"]),
                    }
                    results.append(row)
                    p_values.append(row["pvalue"])
                    log.append(
                        f"[{idx+1}/{total}] {g1} vs {g2} — "
                        f"F={row['pseudo_F']:.3f}, p={row['pvalue']:.4f}"
                    )
                except Exception as exc:
                    log.append(f"[{idx+1}/{total}] {g1} vs {g2} — error: {exc}")

                _write_progress(job_id, {
                    "status": "running", "completed": idx + 1, "total": total,
                    "log": log, "results": results,
                })

            # BH correction
            if results and p_values:
                _, q_values, _, _ = multipletests(p_values, method="fdr_bh")
                for row, q in zip(results, q_values):
                    row["qvalue"] = float(q)

            _write_progress(job_id, {
                "status": "complete", "completed": total, "total": total,
                "log": log, "results": results,
            })
        except Exception as exc:
            _write_progress(job_id, {
                "status": "error", "completed": 0, "total": 0,
                "log": [f"Fatal error: {exc}"], "results": [],
            })

    t = threading.Thread(target=_worker, daemon=True)
    t.start()


def read_permanova_progress(job_id: str) -> dict | None:
    """Read the progress JSON for a background pairwise PERMANOVA job."""
    path = _permanova_progress_path(job_id)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def cleanup_permanova_progress(job_id: str):
    """Remove the progress file for a completed job."""
    path = _permanova_progress_path(job_id)
    path.unlink(missing_ok=True)
