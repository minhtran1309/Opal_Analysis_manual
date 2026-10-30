---
name: analyzing-qupath-cells
description: Use when analysing QuPath cell-detection / cell-measurement exports (Opal / Polaris multiplex IF, TMA or whole-slide) in this project - clustering cells by marker intensity, building an AnnData object, UMAP/Leiden, or mapping clusters back onto tissue coordinates, as in Skin_1_analysis.ipynb.
---

# Analysing QuPath cell exports

## Overview
`scripts/run_pipeline.py` is the `Skin_1_analysis.ipynb` pipeline as a reusable script. It works for any QuPath export that has `Centroid X µm` / `Centroid Y µm` columns and per-marker `Cell: <marker> (...) mean` columns.

Steps, in notebook order, with QC added (see Quality control):
1. Load the export(s), then run the QC filters.
2. Spatial maps of log(mean + 1) intensity for each marker.
3. Keep the whole-cell mean intensities (autofluorescence dropped) and the morphology columns.
4. Build an AnnData object, save the raw matrix, then clip each marker at a high quantile.
5. Boxplot and highest-expressed plots, then `log1p`.
6. Drop the nuclear marker, then PCA → neighbours → UMAP → Leiden.
7. Dotplot per cluster, and cluster maps on the tissue coordinates.

## Running it
Use the project's uv `.venv` (Python 3.12). The notebook's old conda env isn't needed.

```bash
# from the project root
.venv/bin/python .claude/skills/analyzing-qupath-cells/scripts/run_pipeline.py \
    "Cell_dectection_results/*_cell_detection_TMA*.txt" \
    --out results/tma_pipeline --group-by Parent
```

Several files are combined into one object, and `obs['sample']` holds each file's stem. Run with `--help` to see all options.

| Option | Default | When to change |
|---|---|---|
| `--markers` | every `Cell: … mean` column | to analyse a subset of markers |
| `--exclude` | `AF Autofluorescence` | if other channels should be dropped |
| `--nuclear-marker` | `DAPI` | if the panel uses a different nuclear stain |
| `--group-by` | none (one map per file) | `Parent` facets the cluster maps by TMA core |
| `--resolution` | `0.07` (notebook value) | gives coarse clusters; raise it (e.g. 0.3–1.0) for cell types |
| `--scale` | off (notebook didn't scale) | z-scores markers so bright channels don't dominate PCA |
| `--n-pcs` / `--n-neighbors` | 5 / 15 | n-pcs is capped at n_markers − 1 |

## Quality control
QC draws on `tma_exploration.ipynb` and the IMC example `0921_scRNA_IMC.ipynb`. Quantile thresholds are worked out separately for each sample (file), because slides differ. `--no-qc` skips the filters and reproduces the Skin_1 notebook exactly.

| Check | Rule (default) | Option | Figure / table |
|---|---|---|---|
| Missing values | NaN in a marker or a coordinate | always on | `qc_summary.csv` |
| Cell size | keep `Cell: Area` within q0.01–q0.99; small = debris/fragments, big = merged cells | `--area-quantiles` (IMC used 0.05 0.95) | `qc_cell_area.png`, `qc_metrics.png` |
| Low nuclear stain | DAPI below q0.01 | `--min-nuclear-quantile` | `qc_metrics.png` |
| No signal across the panel | fewer than 1 non-DAPI marker above its q0.10 background | `--min-markers` (IMC used 2), `--signal-quantile` | `qc_markers_per_cell.png` |
| High autofluorescence | above the given quantile; off by default | `--max-af-quantile 0.99` | `qc_metrics.png` |
| Small cores | fewer passing cells than N; off by default | `--min-cells-per-group` (needs `--group-by`) | `qc_cells_per_group.png` |
| Where removed cells are | review by eye: clusters of removed cells at edges or folds mean a tissue artefact | — | `qc_removed_<sample>.png` |
| Very bright cells | clip each marker at q0.99 (lower bound 0) before `log1p` | `--clip-quantile` (IMC used 0.95; 0 = off), `--clip-marker M=Q` per marker | `qc_clip_candidates.png/.csv`, `qc_clipping.png` (before/after), `qc_clip_values.csv` |
| Staining batches | per-core median of each raw marker | — | `qc_marker_medians.png/.csv` |
| Batch-driven clusters | each cluster's make-up by core/sample | — | `qc_cluster_composition.png/.csv` |

Why "no signal" uses a background quantile rather than `> 0`: Opal mean intensities are almost never exactly 0. So the IMC rule `sc.pp.filter_cells(min_genes=2)` would remove almost nothing here.

Why clipping matters: a handful of very bright cells (saturation, edge artefacts, over-staining) stretch each marker's range and end up driving PCA and clustering.

### Choosing clip levels for each marker
1. Run with `--qc-only`. It stops after QC and clipping (~15 s on the 3 TMAs; the full run takes minutes):
   ```bash
   .venv/bin/python .claude/skills/analyzing-qupath-cells/scripts/run_pipeline.py \
       "Cell_dectection_results/*_cell_detection_TMA*.txt" --out results/tma_qc --group-by Parent --qc-only
   ```
2. Open `figures/qc_clip_candidates.png`. It shows each marker's histogram (log counts) with dashed lines at the `--clip-candidates` quantiles (default 0.95 0.99 0.995 0.999) and the chosen cap as a solid red line. `qc_clip_candidates.csv` gives the intensity at each candidate and the maximum.
3. Decide for each marker:
   - **Thin, smoothly falling tail** (e.g. p62, MCT4, HCAR1): clip where the tail thins out. q0.99 is usually enough.
   - **Distinct bright population** (a second peak or flat shelf, e.g. LC3AB, Hpyl, the CD147 shoulder): the cap should sit *above* that population, or the marker shouldn't be clipped at all (`M=1`). Otherwise the population collapses into one value. Check in QuPath whether it's real staining or an artefact before keeping it.
   - **A few isolated extreme cells** far past everything else: clip below them.
4. Re-run with the overrides until `qc_clipping.png` looks right. Then drop `--qc-only`:
   ```bash
   ... --clip-quantile 0.99 --clip-marker Hpyl=0.9995 LC3AB=1 CD147=0.995
   ```
The chosen levels are stored in `uns['qc_params']['clip_marker']` and the cap values in `uns['clip_values']`.

## Outputs (in `--out`)
- `qc_summary.csv`: for each sample/core, the cell count, how many cells each check removed, how many passed, and the % removed. Also `qc_clip_values.csv`, `qc_marker_medians.csv` and `qc_cluster_composition.csv`.
- `cells_raw.h5ad`: cells that passed QC, with unclipped intensities in `X`, morphology and IDs in `obs`, coordinates in `obsm['spatial']`, and the QC settings in `uns['qc_params']`.
- `cells_clustered.h5ad`: `log1p` of the clipped values in `X`, raw values in `layers['raw']`, clipped values in `layers['clipped']`, plus `obs['leiden']`, `obsm['X_pca']`, `obsm['X_umap']`, and the run parameters in `uns['leiden_params']`.
- `figures/`: the `qc_*.png` files listed under Quality control, plus `spatial_<marker>.png`, `marker_boxplot.png`, `highest_expr_*.png`, `pca.png`, `umap_leiden.png`, `dotplot_leiden.png`, `spatial_leiden_<sample>.png`.

## Using the results in a notebook
```python
import anndata as ad
adata = ad.read_h5ad("results/tma_pipeline/cells_clustered.h5ad")
```
After that, continue with custom analysis (e.g. squidpy neighbourhood enrichment, using `library_key="Parent"` so cores aren't linked to each other).

## Common mistakes
- **Leaving the resolution at 0.07 and reading the clusters as cell types.** The default reproduces the example; tune it and check the dotplot.
- **Comparing samples without `--scale` or batch correction.** Differences in intensity between slides end up driving the clusters. Check `umap_leiden.png` coloured by sample.
- **Marker names come from the text before ` (`.** So `Cell: DAPI (DAPI) mean` and `Cell: DAPI mean` both become `DAPI`. Pass `--markers` using these short names.
- **Comma-separated exports** (e.g. `TMA2_by_parent/*.csv`) are detected by the `.csv` extension. Don't mix them with the full `.txt` export of the same TMA, or those cells are counted twice.
