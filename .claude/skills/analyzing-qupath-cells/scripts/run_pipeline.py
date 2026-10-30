"""QuPath cell-measurement pipeline, generalised from Skin_1_analysis.ipynb.

Steps (notebook order, with QC added):
  1. load QuPath export(s)
  QC. flag/remove bad cells (missing values, cell-area outliers, low nuclear stain,
      no signal across the panel, optional high autofluorescence / small groups),
      write QC report + figures; later clip each marker at a high quantile
  2. spatial maps of log(mean+1) intensity per marker
  3. select whole-cell mean intensities (minus autofluorescence) + morphology
  4. build AnnData, save raw matrix; per-group marker medians (batch check)
  5. marker boxplot, highest-expressed plot, log1p
  6. drop nuclear stain -> PCA -> neighbours -> UMAP -> Leiden
  7. dotplot per cluster, cluster maps on tissue, cluster composition per group
"""

import argparse
import glob
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

X_COL, Y_COL = "Centroid X µm", "Centroid Y µm"
AREA_COL = "Cell: Area"
# "Cell: CD8 (Opal 480) mean" or "Cell: DAPI mean" -> marker name
MEAN_RE = re.compile(r"^Cell: (?P<marker>[^(]+?)(?: \(.*\))? mean$")
SHAPE_KEYS = ("Area", "Perimeter", "Circularity", "caliper", "Eccentricity", "area ratio")
AF_NAMES = ("AF", "Autofluorescence")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+", help="QuPath export files or glob patterns (.txt/.tsv tab-separated, .csv)")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--markers", nargs="+", help="markers to keep (default: every 'Cell: ... mean' column)")
    p.add_argument("--exclude", nargs="+", default=list(AF_NAMES),
                   help="markers to drop entirely (default: autofluorescence)")
    p.add_argument("--nuclear-marker", default="DAPI", help="kept in the matrix but left out of clustering")
    p.add_argument("--spatial-markers", nargs="+", help="markers for the step-2 intensity maps (default: all)")
    p.add_argument("--group-by", help="obs column to facet plots/QC by, e.g. Parent for TMA cores")
    p.add_argument("--scale", action="store_true", help="z-score markers before PCA (not in the original notebook)")
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--n-pcs", type=int, default=5)
    p.add_argument("--resolution", type=float, default=0.07, help="Leiden resolution (notebook used 0.07)")
    qc = p.add_argument_group("quality control (quantiles are computed per sample)")
    qc.add_argument("--no-qc", action="store_true", help="skip filtering, reproduce the notebook exactly")
    qc.add_argument("--area-quantiles", nargs=2, type=float, default=[0.01, 0.99], metavar=("LO", "HI"),
                    help="keep cells with Cell: Area inside these quantiles")
    qc.add_argument("--min-nuclear-quantile", type=float, default=0.01,
                    help="drop cells with nuclear-marker mean below this quantile")
    qc.add_argument("--max-af-quantile", type=float,
                    help="drop cells with autofluorescence above this quantile, e.g. 0.99 (off by default)")
    qc.add_argument("--signal-quantile", type=float, default=0.10,
                    help="a marker counts as expressed when above this per-sample quantile (background)")
    qc.add_argument("--min-markers", type=int, default=1,
                    help="drop cells expressing fewer panel markers than this (IMC example used 2)")
    qc.add_argument("--min-cells-per-group", type=int, default=0,
                    help="drop --group-by groups (cores) with fewer passing cells than this")
    qc.add_argument("--clip-quantile", type=float, default=0.99,
                    help="clip each marker at this quantile before log1p (IMC example used 0.95; 0 = off)")
    qc.add_argument("--clip-marker", nargs="+", default=[], metavar="MARKER=Q",
                    help="per-marker clip quantile overriding --clip-quantile, e.g. Hpyl=0.999 LC3AB=1 (1 = no clip)")
    qc.add_argument("--clip-candidates", nargs="+", type=float, default=[0.95, 0.99, 0.995, 0.999],
                    help="quantiles drawn on qc_clip_candidates.png to help choose the clip level")
    qc.add_argument("--qc-only", action="store_true",
                    help="stop after QC and clipping figures (fast loop for choosing thresholds)")
    args = p.parse_args()
    try:
        args.clip_marker = {m: float(q) for m, q in (s.split("=", 1) for s in args.clip_marker)}
    except ValueError:
        p.error("--clip-marker expects MARKER=QUANTILE, e.g. Hpyl=0.999")
    return args


def load(inputs):
    files = sorted({Path(f) for pat in inputs for f in (glob.glob(pat) or [pat])})
    frames = []
    for f in files:
        sep = "," if f.suffix.lower() == ".csv" else "\t"
        df = pd.read_csv(f, sep=sep, low_memory=False)
        frames.append(df.assign(sample=f.stem))
        print(f"loaded {f} ({len(df)} cells)")
    return pd.concat(frames, ignore_index=True)


def mean_columns(raw):
    return {m.group("marker"): c for c in raw.columns if (m := MEAN_RE.match(c))}


def marker_columns(found, markers, exclude):
    if markers:
        missing = [m for m in markers if m not in found]
        if missing:
            raise SystemExit(f"markers not found: {missing}; available: {list(found)}")
        found = {m: found[m] for m in markers}
    return {m: c for m, c in found.items() if m not in exclude}


def scatter(ax, xy, c, title, **kw):
    s = ax.scatter(xy[:, 0], xy[:, 1], c=c, s=1, rasterized=True, **kw)
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.invert_yaxis()  # image coordinates: y increases downwards
    return s


def facet_grid(n, panel=4):
    ncols = min(n, 6)
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel * ncols, panel * nrows), squeeze=False)
    for ax in axes.flat[n:]:
        ax.axis("off")
    return fig, axes


def save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


def quality_control(raw, all_means, cols, args, out):
    """Flag failing cells, write qc_summary.csv + figures, return the passing cells."""
    figs = out / "figures"
    per_sample = raw.groupby("sample")

    def outside(col, lo=None, hi=None):
        v = raw[col]
        bad = pd.Series(False, index=raw.index)
        if lo is not None:
            bad |= v < per_sample[col].transform(lambda s: s.quantile(lo))
        if hi is not None:
            bad |= v > per_sample[col].transform(lambda s: s.quantile(hi))
        return bad

    # thresholds: (column, lo quantile, hi quantile) for plotting
    checks = {AREA_COL: tuple(args.area_quantiles)}
    fails = pd.DataFrame(index=raw.index)
    fails["missing_values"] = raw[list(cols.values()) + [X_COL, Y_COL]].isna().any(axis=1)
    fails["cell_area"] = outside(AREA_COL, *args.area_quantiles)
    if args.nuclear_marker in all_means:
        col = all_means[args.nuclear_marker]
        fails["low_nuclear"] = outside(col, lo=args.min_nuclear_quantile)
        checks[col] = (args.min_nuclear_quantile, None)
    af_col = next((all_means[n] for n in AF_NAMES if n in all_means), None)
    if af_col and args.max_af_quantile:
        fails["high_autofluorescence"] = outside(af_col, hi=args.max_af_quantile)
        checks[af_col] = (None, args.max_af_quantile)
    # cells with no signal across the panel: Opal means are rarely exactly 0, so "expressed"
    # means above the per-sample background quantile of that marker
    panel = [c for m, c in cols.items() if m != args.nuclear_marker]
    background = per_sample[panel].transform("quantile", args.signal_quantile)
    n_expressed = (raw[panel] > background).sum(axis=1)
    fails["no_signal"] = n_expressed < args.min_markers

    facet = args.group_by or "sample"
    keys = ["sample", facet] if args.group_by else ["sample"]
    passing = ~fails.any(axis=1)
    if args.min_cells_per_group and args.group_by:
        n_pass = passing.groupby([raw[k] for k in keys]).transform("sum")
        fails["small_group"] = passing & (n_pass < args.min_cells_per_group)
        passing &= ~fails["small_group"]

    summary = fails.groupby([raw[k] for k in keys]).sum()
    summary.insert(0, "n_cells", raw.groupby(keys).size())
    summary["n_pass"] = passing.groupby([raw[k] for k in keys]).sum()
    summary["pct_removed"] = (100 * (1 - summary["n_pass"] / summary["n_cells"])).round(1)
    summary.to_csv(out / "qc_summary.csv")
    print(f"wrote {out / 'qc_summary.csv'}")
    total = summary.drop(columns="pct_removed").sum()
    print(f"QC: {int(total['n_pass'])}/{int(total['n_cells'])} cells pass; "
          + ", ".join(f"{k}={int(total[k])}" for k in fails.columns))

    # metric distributions per sample with the per-sample thresholds
    metrics = [AREA_COL, "Nucleus: Area", "Nucleus/Cell area ratio"] + [c for c in checks if c != AREA_COL]
    metrics = [m for m in dict.fromkeys(metrics) if m in raw.columns]
    if af_col and af_col not in metrics:
        metrics.append(af_col)
    fig, axes = plt.subplots(len(metrics), 1, figsize=(max(6, 2 * raw["sample"].nunique()), 3 * len(metrics)),
                             squeeze=False)
    for ax, m in zip(axes.flat, metrics):
        sns.violinplot(data=raw, x="sample", y=m, cut=0, inner="quartile", ax=ax)
        for i, (s, grp) in enumerate(raw.groupby("sample", sort=False)):
            for q in filter(lambda q: q is not None, checks.get(m, ())):
                ax.hlines(grp[m].quantile(q), i - 0.4, i + 0.4, colors="red", linestyles="--")
        ax.set_xlabel("")
    fig.suptitle("QC metrics per sample (red dashes = filter thresholds)")
    save(fig, figs / "qc_metrics.png")

    # cell size distribution: too small = fragments/debris, too big = merged cells
    fig, ax = plt.subplots(figsize=(8, 4))
    for i, (s, grp) in enumerate(raw.groupby("sample", sort=False)):
        color = f"C{i}"
        ax.hist(grp[AREA_COL], bins=np.geomspace(raw[AREA_COL].clip(lower=1e-2).min(), raw[AREA_COL].max(), 100),
                histtype="step", color=color, label=s)
        for q in args.area_quantiles:
            ax.axvline(grp[AREA_COL].quantile(q), color=color, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel(f"{AREA_COL} (µm²)")
    ax.set_ylabel("cells")
    ax.set_title("Cell size (dashes = per-sample area cut-offs)")
    ax.legend(fontsize=7)
    save(fig, figs / "qc_cell_area.png")

    # how many panel markers each cell expresses above background
    fig, ax = plt.subplots(figsize=(8, 4))
    pd.crosstab(n_expressed, raw["sample"], normalize="columns").plot.bar(ax=ax)
    ax.axvline(args.min_markers - 0.5, color="red", ls="--", label="min markers")
    ax.set_xlabel(f"markers above the per-sample {args.signal_quantile:.0%} quantile")
    ax.set_ylabel("fraction of cells")
    ax.legend(fontsize=7)
    save(fig, figs / "qc_markers_per_cell.png")

    # cells per group, before and after QC
    counts = summary[["n_cells", "n_pass"]].reset_index()
    counts["label"] = counts[keys].astype(str).agg(" | ".join, axis=1) if len(keys) > 1 else counts["sample"]
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(counts)), 4))
    ax.bar(counts["label"], counts["n_cells"], color="lightgrey", label="all cells")
    ax.bar(counts["label"], counts["n_pass"], color="tab:blue", label="pass QC")
    if args.min_cells_per_group:
        ax.axhline(args.min_cells_per_group, color="red", ls="--", label="min cells per group")
    ax.tick_params(axis="x", rotation=90)
    ax.set_ylabel("cells")
    ax.legend()
    save(fig, figs / "qc_cells_per_group.png")

    # where the removed cells sit: edges, folds or whole cores are artefacts
    for s in raw["sample"].unique():
        sub = raw[raw["sample"] == s]
        groups = sorted(sub[facet].astype(str).unique())
        fig, axes = facet_grid(len(groups)) if args.group_by else plt.subplots(figsize=(10, 10), squeeze=False)
        for ax, g in zip(axes.flat, groups):
            gsub = sub[sub[facet].astype(str) == g]
            ok = passing[gsub.index]
            xy = gsub[[X_COL, Y_COL]].to_numpy()
            scatter(ax, xy[ok.values], "lightgrey", f"{g} ({(~ok).sum()} removed)")
            ax.scatter(xy[~ok.values, 0], xy[~ok.values, 1], c="red", s=2, rasterized=True)
        save(fig, figs / f"qc_removed_{s}.png")

    return raw[passing].reset_index(drop=True)


def main():
    args = parse_args()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    raw = load(args.inputs)
    all_means = mean_columns(raw)
    cols = marker_columns(all_means, args.markers, args.exclude)
    markers = list(cols)
    print(f"markers: {markers}")
    unknown = set(args.clip_marker) - set(markers)
    if unknown:
        raise SystemExit(f"--clip-marker: unknown markers {sorted(unknown)}; available: {markers}")
    if args.group_by:
        raw[args.group_by] = raw[args.group_by].astype(str)
    if not args.no_qc:
        raw = quality_control(raw, all_means, cols, args, out)

    # 2. spatial intensity maps, one figure per marker with a panel per sample
    samples = raw["sample"].unique()
    for m in args.spatial_markers or markers:
        fig, axes = plt.subplots(1, len(samples), figsize=(8 * len(samples), 8), squeeze=False)
        for ax, s in zip(axes.flat, samples):
            sub = raw[raw["sample"] == s]
            sc_ = scatter(ax, sub[[X_COL, Y_COL]].to_numpy(), np.log1p(sub[cols[m]]), s, cmap="RdBu")
            fig.colorbar(sc_, ax=ax, shrink=0.6, label=f"log({m} mean + 1)")
        save(fig, out / "figures" / f"spatial_{m}.png")

    # 3-4. morphology + identifiers in obs, marker intensities in X
    shape_cols = [c for c in raw.columns if "µm" in c or any(k in c for k in SHAPE_KEYS)]
    id_cols = [c for c in ("sample", "Image", "Parent", "Classification", "Class") if c in raw.columns]
    obs = raw[id_cols + shape_cols].copy()
    for c in id_cols:
        obs[c] = obs[c].astype(str).astype("category")
    obs.columns = obs.columns.str.replace("/", "-")  # h5ad keys can't contain '/'
    obs.index = (raw["sample"] + "_" + raw.index.astype(str)).values

    adata = ad.AnnData(X=raw[list(cols.values())].to_numpy(np.float32), obs=obs,
                       var=pd.DataFrame(index=markers))
    adata.obsm["spatial"] = raw[[X_COL, Y_COL]].to_numpy()
    adata.layers["raw"] = adata.X.copy()
    adata.uns["qc_params"] = {"applied": not args.no_qc, "area_quantiles": args.area_quantiles,
                              "min_nuclear_quantile": args.min_nuclear_quantile,
                              "max_af_quantile": args.max_af_quantile or 0,
                              "signal_quantile": args.signal_quantile, "min_markers": args.min_markers,
                              "min_cells_per_group": args.min_cells_per_group,
                              "clip_quantile": args.clip_quantile, "clip_marker": args.clip_marker}
    adata.write_h5ad(out / "cells_raw.h5ad", compression="gzip")
    print(f"wrote {out / 'cells_raw.h5ad'} {adata.shape}")

    # clip very bright cells so a few saturated/over-stained cells don't drive PCA and clustering
    df = adata.to_df()
    clip_q = pd.Series({m: args.clip_marker.get(m, args.clip_quantile) for m in markers})

    # candidate levels: histogram per marker with each candidate quantile drawn on it, so a thin
    # tail (clip it) can be told apart from a distinct bright population (keep it: raise the quantile)
    cands = sorted(set(args.clip_candidates))
    table = pd.DataFrame({f"q{q}": df.quantile(q) for q in cands})
    table["max"] = df.max()
    table.insert(0, "chosen_quantile", clip_q)
    table.to_csv(out / "qc_clip_candidates.csv")
    ncols = min(len(markers), 4)
    nrows = -(-len(markers) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows), squeeze=False)
    colors = plt.get_cmap("viridis")(np.linspace(0, 0.9, len(cands)))
    for ax, m in zip(axes.flat, markers):
        ax.hist(df[m], bins=100, log=True, color="lightgrey")
        for q, c in zip(cands, colors):
            ax.axvline(df[m].quantile(q), color=c, ls="--", lw=1, label=f"q{q} = {df[m].quantile(q):.2f}")
        if 0 < clip_q[m] < 1:
            ax.axvline(df[m].quantile(clip_q[m]), color="red", lw=2, label=f"chosen q{clip_q[m]}")
        ax.set_title(f"{m} (max {df[m].max():.2f})")
        ax.legend(fontsize=7, title="quantile = intensity", title_fontsize=7)
    for ax in axes.flat[len(markers):]:
        ax.axis("off")
    fig.suptitle("Choosing clip levels: clip thin tails, keep distinct bright populations")
    save(fig, out / "figures" / "qc_clip_candidates.png")

    caps = pd.Series({m: df[m].quantile(q) if 0 < q < 1 else np.inf for m, q in clip_q.items()})
    if np.isfinite(caps).any():
        clipped = df.clip(lower=0, upper=caps, axis=1)
        pd.DataFrame({"quantile": clip_q, "clip_value": caps, "max_before": df.max(),
                      "n_clipped": (df > caps).sum()}).to_csv(out / "qc_clip_values.csv")
        fig, axes = plt.subplots(2, len(markers), figsize=(3 * len(markers), 6), squeeze=False)
        for j, m in enumerate(markers):
            axes[0, j].hist(df[m], bins=50, log=True)
            if np.isfinite(caps[m]):
                axes[0, j].axvline(caps[m], color="red", ls="--")
            axes[0, j].set_title(f"{m} (q{clip_q[m]})" if np.isfinite(caps[m]) else f"{m} (not clipped)")
            axes[1, j].hist(clipped[m], bins=50)
        axes[0, 0].set_ylabel("before clipping (log count)")
        axes[1, 0].set_ylabel("after clipping")
        save(fig, out / "figures" / "qc_clipping.png")
        adata.X = clipped.to_numpy(np.float32)
        adata.layers["clipped"] = adata.X.copy()
        adata.uns["clip_values"] = caps[np.isfinite(caps)].to_dict()

    if args.qc_only:
        print("--qc-only: stopping after QC; inspect figures/qc_*.png and qc_*.csv")
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    adata.to_df().boxplot(ax=ax, rot=90)
    save(fig, out / "figures" / "marker_boxplot.png")

    # batch check: per-group median intensity, scaled per marker; a group that is bright
    # across all markers points to staining/exposure differences rather than biology
    facet = args.group_by or "sample"
    medians = adata.to_df(layer="raw").groupby(adata.obs[facet].values).median()
    medians.to_csv(out / "qc_marker_medians.csv")
    scaled = (medians - medians.min()) / (medians.max() - medians.min())
    fig, ax = plt.subplots(figsize=(0.6 * len(markers) + 2, 0.3 * len(medians) + 2))
    sns.heatmap(scaled, cmap="viridis", ax=ax, cbar_kws={"label": "median, scaled per marker"})
    ax.set_ylabel(facet)
    save(fig, out / "figures" / "qc_marker_medians.png")

    # 5. log1p
    ax = sc.pl.highest_expr_genes(adata, n_top=len(markers), show=False)
    save(ax.figure, out / "figures" / "highest_expr_raw.png")
    sc.pp.log1p(adata)
    ax = sc.pl.highest_expr_genes(adata, n_top=len(markers), show=False)
    save(ax.figure, out / "figures" / "highest_expr_log1p.png")

    # 6. cluster on the non-nuclear markers
    clust = adata[:, [m for m in markers if m != args.nuclear_marker]].copy()
    if args.scale:
        sc.pp.scale(clust, max_value=10)
    n_pcs = min(args.n_pcs, clust.n_vars - 1)
    sc.tl.pca(clust, n_comps=n_pcs, svd_solver="arpack")
    sc.pl.pca(clust, show=False)
    save(plt.gcf(), out / "figures" / "pca.png")
    sc.pp.neighbors(clust, n_neighbors=args.n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(clust)
    sc.tl.leiden(clust, resolution=args.resolution, flavor="igraph", n_iterations=2)
    print(clust.obs["leiden"].value_counts().sort_index().to_string())

    sc.pl.umap(clust, color=["leiden", "sample"], size=8, show=False)
    save(plt.gcf(), out / "figures" / "umap_leiden.png")
    sc.pl.dotplot(clust, clust.var_names, groupby="leiden", show=False)
    save(plt.gcf(), out / "figures" / "dotplot_leiden.png")

    adata.obs["leiden"] = clust.obs["leiden"]
    adata.obsm["X_pca"] = clust.obsm["X_pca"]
    adata.obsm["X_umap"] = clust.obsm["X_umap"]
    adata.uns["leiden_params"] = {"resolution": args.resolution, "n_neighbors": args.n_neighbors,
                                  "n_pcs": n_pcs, "scaled": args.scale}

    # 7. clusters back on tissue
    cats = adata.obs["leiden"].cat.categories
    cmap = dict(zip(cats, map(matplotlib.colors.to_hex, plt.get_cmap("tab20").colors)))
    for s in samples:
        sub = adata[adata.obs["sample"] == s]
        groups = sorted(sub.obs[facet].unique())
        fig, axes = facet_grid(len(groups)) if args.group_by else plt.subplots(figsize=(10, 10), squeeze=False)
        for ax, g in zip(axes.flat, groups):
            gsub = sub[sub.obs[facet] == g]
            scatter(ax, gsub.obsm["spatial"], gsub.obs["leiden"].map(cmap).tolist(), f"{g} (n={gsub.n_obs})")
        fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=cmap[c], label=c) for c in cats],
                   title="leiden", loc="center right")
        save(fig, out / "figures" / f"spatial_leiden_{s}.png")

    # cluster composition per group: a cluster drawn almost entirely from one core/sample
    # is more likely a batch effect than a cell type
    comp = pd.crosstab(adata.obs[facet], adata.obs["leiden"], normalize="index")
    comp.to_csv(out / "qc_cluster_composition.csv")
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(comp)), 4))
    comp.plot.bar(stacked=True, ax=ax, color=[cmap[c] for c in comp.columns], width=0.9)
    ax.set_ylabel("fraction of cells")
    ax.legend(title="leiden", bbox_to_anchor=(1, 1))
    save(fig, out / "figures" / "qc_cluster_composition.png")

    adata.write_h5ad(out / "cells_clustered.h5ad", compression="gzip")
    print(f"wrote {out / 'cells_clustered.h5ad'}")


if __name__ == "__main__":
    main()
