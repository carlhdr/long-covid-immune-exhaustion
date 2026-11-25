"""
PCA of the GSE275334 immune-exhaustion panel, with age and sex as candidate
confounders.

Reads data/samples.csv and data/expr_log2.csv (written by load_data.py).
Writes:
  data/pca_scores.csv    per-sample PC scores and the metadata used
  figures/pca_groups.png the figure

Method: genes are z-scored across samples, so each of the 635 genes contributes
equally; PCA is then an SVD of the centred sample x gene matrix. Group, age and
sex associations with each of the first three PCs are tested non-parametrically
(Kruskal-Wallis across the three groups, Spearman for age, Mann-Whitney for
sex). These p-values are uncorrected and exploratory.

Run:  python src/pca.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))  # so style.py imports
from style import apply_figure_style, panel_letter, strip_with_median  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIGS = ROOT / "figures"

GROUPS = ["Healthy control", "Long COVID", "ME/CFS"]
COLORS = {
    "Healthy control": "#8c8c8c",   # neutral grey: the reference group
    "Long COVID": "#1f78b4",        # blue
    "ME/CFS": "#ff7f00",            # orange; blue/orange stays CVD-safe (S4.5)
}


def run_pca(expr: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Z-score genes, then SVD. Returns (scores [samples x PC], variance ratio)."""
    z = ((expr.T - expr.mean(axis=1)) / expr.std(axis=1, ddof=1)).values
    u, s, _ = np.linalg.svd(z - z.mean(axis=0), full_matrices=False)
    return u * s, s**2 / np.sum(s**2)


def associations(scores: np.ndarray, samples: pd.DataFrame, n_pcs: int = 3) -> pd.DataFrame:
    rows = []
    for i in range(n_pcs):
        pc = scores[:, i]
        rows.append(
            {
                "pc": f"PC{i + 1}",
                "group_kruskal_p": stats.kruskal(
                    *[pc[samples.group.values == g] for g in GROUPS]
                ).pvalue,
                "age_spearman_rho": stats.spearmanr(pc, samples.age).statistic,
                "age_spearman_p": stats.spearmanr(pc, samples.age).pvalue,
                "sex_mannwhitney_p": stats.mannwhitneyu(
                    pc[samples.sex.values == "Female"], pc[samples.sex.values == "Male"]
                ).pvalue,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    samples = pd.read_csv(DATA / "samples.csv")
    expr = pd.read_csv(DATA / "expr_log2.csv", index_col=0)[samples.sample_id]
    assert list(expr.columns) == list(samples.sample_id)

    scores, var = run_pca(expr)
    flagged = samples.qc_issues.notna().values

    out = samples.copy()
    for i in range(5):
        out[f"PC{i + 1}"] = scores[:, i]
    out.to_csv(DATA / "pca_scores.csv", index=False)

    assoc = associations(scores, samples)
    print(f"variance explained, PC1-5 (%): {np.round(var[:5] * 100, 1)}")
    print(assoc.round(3).to_string(index=False))

    # ---- figure ---------------------------------------------------------
    apply_figure_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2), width_ratios=[1.35, 1])

    # (a) PC1 vs PC2. QC-flagged samples drawn as open markers (S1.1): they are
    # kept in the PCA but must be visually distinguishable.
    for g in GROUPS:
        m = (samples.group.values == g) & ~flagged
        ax1.scatter(scores[m, 0], scores[m, 1], s=34, c=COLORS[g], label=g,
                    edgecolors="white", linewidths=0.5, zorder=3)
    m = flagged
    ax1.scatter(scores[m, 0], scores[m, 1], s=40, facecolors="none",
                edgecolors=[COLORS[g] for g in samples.group.values[m]],
                linewidths=1.3, label="flagged by vendor QC", zorder=4)

    ax1.set_xlabel(f"PC1 ({var[0] * 100:.0f}% of variance)")
    ax1.set_ylabel(f"PC2 ({var[1] * 100:.0f}% of variance)")
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.margins(0.08)
    ax1.legend(loc="upper left", bbox_to_anchor=(0.0, 1.02), handletextpad=0.3,
               borderpad=0.0, labelspacing=0.3)
    panel_letter(ax1, "a", dx=-0.10)

    # (b) PC2 is the only one of the first three with any group association.
    strip_with_median(
        ax2,
        groups=samples.group.values,
        values=scores[:, 1],
        colors=COLORS,
    )
    p2 = assoc.loc[assoc.pc == "PC2", "group_kruskal_p"].item()
    ax2.set_ylabel("PC2 score")
    ax2.set_xticks(range(len(GROUPS)))
    ax2.set_xticklabels(["Healthy\ncontrol", "Long\nCOVID", "ME/CFS"])
    ax2.annotate(f"Kruskal-Wallis $p$ = {p2:.2f}\n(uncorrected, 3 PCs tested)",
                 xy=(0.5, 0.01), xycoords="axes fraction", ha="center",
                 va="bottom", fontsize=plt.rcParams["legend.fontsize"],
                 color="#4d4d4d")
    ax2.margins(y=0.12)
    panel_letter(ax2, "b", dx=-0.22)

    fig.suptitle(
        "A 635-gene immune-exhaustion panel does not separate Long COVID,\n"
        "ME/CFS and healthy controls in PBMCs",
        y=1.06, ha="center",
    )
    FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / "pca_groups.png", bbox_inches="tight")
    print(f"wrote {FIGS / 'pca_groups.png'}")


if __name__ == "__main__":
    main()
