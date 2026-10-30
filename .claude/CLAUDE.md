# Conver_cell_detection

Analysis of QuPath cell-detection exports from Opal multiplex immunofluorescence (Polaris) images, done with scanpy / squidpy / AnnData.

## Environment
- uv-managed `.venv`, Python 3.12 (`pyproject.toml`, `uv.lock`). Run code with `.venv/bin/python` or `uv run`, and add packages with `uv add <pkg>`.
- Notebooks must use the registered kernel **Python (conver-cell-detection)** (kernelspec name `conver-cell-detection`), which points at `.venv/bin/python3` by absolute path. Never save or generate notebooks with the generic `python3` kernel: it runs whichever `python` is first on PATH (Anaconda's on this machine) and fails with `No module named 'anndata'`. When running notebooks with nbconvert, pass `--ExecutePreprocessor.kernel_name=conver-cell-detection`. Setup and troubleshooting are in `README.md`.
- `Skin_1_analysis.ipynb` came from an older Python 3.8 conda env (`imctools2`) and reads from `/Volumes/BiomedML/...`, which is only available when that share is mounted. Treat it as a reference, not something to run here.
- This folder isn't a git repository.

## Data
- `Cell_dectection_results/YYYYMMDD_cell_detection_TMA{1,2,3}.txt`: tab-separated QuPath exports, ~80–100k cells per TMA.
  - `Parent` holds the core ID (e.g. `COL-502BG`).
  - Coordinates are in `Centroid X µm` / `Centroid Y µm`.
  - Each marker has mean-intensity columns named `{Nucleus,Cell,Cytoplasm}: <marker> (Opal NNN) mean`.
  - TMA panel: p62, CD147, LC3AB, Hpyl, MCT4, MCT1, HCAR1, plus DAPI and Autofluorescence. Those last two have no Opal suffix, e.g. `Cell: DAPI mean`.
- `Cell_dectection_results/TMA2_by_parent/*.csv`: TMA2 split per core, comma-separated. These cells duplicate the TMA2 `.txt`, so never load both.
- `results/`: pipeline outputs (`.h5ad`, figures).

## Notebooks
- `Skin_1_analysis.ipynb`: the reference pipeline (one skin section, immune panel CD8 / PD-L1 / FoxP3 / PD-1 / PanCK / CD68).
- `tma_exploration.ipynb`: the skill's pipeline (QC → clipping → clustering → spatial) walked through step by step on one assay (`ASSAY`, a `Parent` value; default `COL-502BG` on TMA1). Saves to `results/assays/<tma>_<assay>.h5ad`. The previous version is in `tma_exploration_backup.ipynb`.
- A TMA holds many small assays (cores), identified by `Parent`. Per-assay cluster labels aren't comparable across assays; use the skill script to cluster assays jointly.
- `data_exploration.py` and `pipeline` are empty placeholders.

## Analysis pipeline
Use the project skill `analyzing-qupath-cells` (`.claude/skills/analyzing-qupath-cells/`). It wraps the Skin_1 pipeline as `scripts/run_pipeline.py`:

```bash
.venv/bin/python .claude/skills/analyzing-qupath-cells/scripts/run_pipeline.py \
    "Cell_dectection_results/*_cell_detection_TMA*.txt" --out results/tma_pipeline --group-by Parent
```

## Conventions
- AnnData layout:
  - `X` holds whole-cell mean intensities, with markers as `var_names` (short names like `CD8`, not full column names).
  - `obs` holds IDs and morphology; `obsm['spatial']` holds centroids in µm; `layers['raw']` keeps the untransformed values.
- Plots on tissue coordinates call `invert_yaxis()`, because image y increases downwards.
- Autofluorescence is always excluded from analysis. DAPI stays in the matrix but is left out of clustering.
- Build spatial graphs per core (`library_key`) so cells in different cores aren't connected.
