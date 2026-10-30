# Conver_cell_detection

Analysis of QuPath cell-detection exports from Opal multiplex immunofluorescence (Polaris) images
with scanpy / squidpy / AnnData.

## Installation (uv)

The project uses [uv](https://docs.astral.sh/uv/) with Python 3.12 (`pyproject.toml`, `uv.lock`, `.python-version`).

### Install uv first

> **Note:** uv is a separate tool and is often **not installed** – it doesn't come with Python,
> Anaconda or macOS, and HPC systems such as Katana don't provide it as a module. Check first:
>
> ```bash
> uv --version        # "command not found" → install it as below
> ```

Install it (no admin rights needed; it goes into `~/.local/bin`):

```bash
# macOS / Linux (including Katana)
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.local/bin/env          # make uv available in this shell (or open a new terminal)
uv --version
```

Alternatives: `brew install uv` (macOS with Homebrew), `pip install --user uv` (anywhere with Python),
or on Windows (PowerShell): `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"`.

If `uv` is still "not found" after installing, `~/.local/bin` isn't on your `PATH`: run
`source ~/.local/bin/env`, or add `export PATH="$HOME/.local/bin:$PATH"` to `~/.bashrc` / `~/.zshrc`.
You don't need to install Python separately: uv downloads Python 3.12 for the project by itself.

### Set up the project

```bash
cd Opal_Analysis_manual          # the cloned project folder

# 1. create .venv and install the locked dependencies (including the dev group: jupyterlab, ipykernel)
uv sync

# 2. register the project kernel with Jupyter  <-- required, see "Jupyter kernel" below
uv run python -m ipykernel install --user \
    --name conver-cell-detection --display-name "Python (conver-cell-detection)"

# 3. check: the kernel should point at this project's .venv
uv run jupyter kernelspec list
```

Step 2 only needs to be repeated if `.venv` is deleted/recreated or the project folder is moved
(the kernel stores the absolute path to `.venv/bin/python3`).

Adding a package: `uv add <package>` (it is available in the notebooks after restarting the kernel).

## Running on Katana (UNSW HPC)

### 1. Log in
```bash
ssh -XYC <zID>@katana2.restech.unsw.edu.au      # e.g. z1234567@katana2.restech.unsw.edu.au
```
`-X`/`-Y` forward X11 so plot windows can open on your machine; `-C` compresses the connection
(helps with X11 over slow links).

### 2. Request an interactive compute session
The login node is shared – don't run analyses on it. Request a compute node with one of these commands:

```bash
# 2 CPUs, 64 GB RAM, 2 h – one assay in the notebook, --qc-only runs
qsub -IX -l select=1:ncpus=2:mem=64gb,walltime=2:00:00

# 2 CPUs, 100 GB RAM, 8 h – full pipeline on all TMAs
qsub -IX -l select=1:ncpus=2:mem=100gb,walltime=8:00:00

# 2 CPUs, 1 GPU, 64 GB RAM, 9 h – GPU work (the current pipeline doesn't need a GPU)
qsub -IX -l select=1:ncpus=2:ngpus=1:mem=64gb,walltime=9:00:00
```

`-I` = interactive, `-X` = X11 forwarding into the job, `select=1` = one node. Wait until the prompt
moves to a compute node, then continue there. `exit` ends the session; `qstat -u $USER` lists your jobs.

### 3. Install uv (once)
uv is **not available on Katana by default** (there is no `uv` module), so `uv` gives
"command not found" until you install it yourself. Check, and install if needed – it goes into your
home directory, no admin rights or modules needed:

```bash
uv --version || curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.local/bin/env          # make uv available now (new logins pick it up automatically)
uv --version
```

If `uv` is still not found in a later session or inside a `qsub` job, add
`export PATH="$HOME/.local/bin:$PATH"` to `~/.bashrc`. See [Install uv first](#install-uv-first) for details.

uv downloads Python 3.12 itself (from `.python-version`), so no Python module has to be loaded.
If your home quota is tight, point uv's cache at your scratch space, e.g. add
`export UV_CACHE_DIR=/srv/scratch/$USER/uv-cache` to `~/.bashrc`.

### 4. Get the project and the data
```bash
git clone https://github.com/minhtran1309/Opal_Analysis_manual.git   # public: no login needed
cd Opal_Analysis_manual
uv sync                                                              # creates .venv
```

The data aren't in the repository. Copy the QuPath exports from your computer (run this **on your
computer**, not on Katana):

```bash
rsync -avP Cell_dectection_results/ <zID>@katana2.restech.unsw.edu.au:Opal_Analysis_manual/Cell_dectection_results/
```

For large transfers, UNSW recommends Katana's data mover node; see the UNSW Research Technology
Services documentation.

### 5. Run the pipeline
Inside the interactive session, from the project folder:

```bash
# fast QC loop: look at results/tma_pipeline/figures/qc_*.png, choose clip levels
uv run python .claude/skills/analyzing-qupath-cells/scripts/run_pipeline.py \
    "Cell_dectection_results/*_cell_detection_TMA*.txt" --out results/tma_pipeline --group-by Parent --qc-only

# full run
uv run python .claude/skills/analyzing-qupath-cells/scripts/run_pipeline.py \
    "Cell_dectection_results/*_cell_detection_TMA*.txt" --out results/tma_pipeline --group-by Parent
```

The script draws no windows: figures go to `results/tma_pipeline/figures/` as PNG files. View them
with an X11 image viewer (e.g. `display file.png`, if available) or copy them back with `rsync`/`scp`.

To run the notebook without a browser (outputs are saved into the notebook file):

```bash
uv run python -m ipykernel install --user --name conver-cell-detection --display-name "Python (conver-cell-detection)"
uv run jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=conver-cell-detection tma_exploration.ipynb
```

## Running notebooks

Open the notebooks with any Jupyter front end (JupyterLab, Anaconda's Jupyter, VS Code) and select the
kernel **Python (conver-cell-detection)**. The notebooks in this project are saved with that kernel.

```bash
uv run jupyter lab          # JupyterLab from the project environment
```

## Jupyter kernel – why it matters, and fixing `No module named 'anndata'`

### Symptom
The first cell of a notebook fails with

```
----> 4 import anndata as ad
ModuleNotFoundError: No module named 'anndata'
```

(or the same error for `scanpy`, `squidpy`, …). In a cell, `import sys; print(sys.executable)` shows
something other than `.../Conver_cell_detection/.venv/bin/python3`, e.g. `/Users/<you>/anaconda3/bin/python`.

### Cause
The notebook is running on the wrong Python. The project's packages are installed only in `.venv`.
A notebook set to the generic **`python3` / "Python 3 (ipykernel)"** kernel does not use `.venv`: that
kernel launches whatever `python` comes first on your `PATH`. On a machine with Anaconda that is
`~/anaconda3/bin/python`, which has none of the project packages. It only "works" when Jupyter itself
happens to be started from inside `.venv`, which is why a notebook can run for one person (or one
launch method) and fail for another.

The registered kernel **`conver-cell-detection`** avoids this: its `kernel.json` contains the
absolute path to `.venv/bin/python3`, so it uses the project environment no matter which Jupyter
(Anaconda, VS Code, `uv run jupyter lab`) opens the notebook.

### Fix
1. Make sure the kernel is registered (installation step 2):
   ```bash
   uv run jupyter kernelspec list
   # conver-cell-detection   /Users/<you>/Library/Jupyter/kernels/conver-cell-detection
   ```
   If it's missing, run the `ipykernel install` command from the installation section.
2. Check it points at this project:
   ```bash
   cat ~/Library/Jupyter/kernels/conver-cell-detection/kernel.json
   # "argv": ["/.../Conver_cell_detection/.venv/bin/python3", "-m", "ipykernel_launcher", ...]
   ```
   If the path is wrong (folder moved, `.venv` recreated), re-run the install command; it overwrites the old spec.
3. Select the kernel in the notebook:
   - **JupyterLab / Jupyter Notebook:** *Kernel → Change Kernel… → Python (conver-cell-detection)*
   - **VS Code:** kernel picker (top right) → *Jupyter Kernel… → Python (conver-cell-detection)*

   Selecting it saves the choice in the notebook's metadata, so it opens with the right kernel next time.
   If a notebook was open while its kernel setting changed on disk, close and reopen it.
4. Verify in the first cell:
   ```python
   import sys; print(sys.executable)   # .../Conver_cell_detection/.venv/bin/python3
   ```

### For notebooks created or generated by scripts
Set the notebook metadata to the project kernel, not `python3`:

```json
"kernelspec": {"display_name": "Python (conver-cell-detection)", "language": "python", "name": "conver-cell-detection"}
```

When executing notebooks from the command line, name the kernel explicitly:

```bash
uv run jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=conver-cell-detection tma_exploration.ipynb
```

## Project layout

| Path | Contents |
|---|---|
| `tma_exploration.ipynb` | step-by-step manual: the pipeline on one TMA assay (QC → clipping → clustering → spatial) |
| `.claude/skills/analyzing-qupath-cells/` | the pipeline as a script (`scripts/run_pipeline.py`) + documentation (`SKILL.md`, incl. all QC steps) |
| `.claude/CLAUDE.md` | project notes for Claude Code |
| `Cell_dectection_results/` | *not in the repository* – put the QuPath exports here, one tab-separated file per TMA (`Parent` = assay/core ID) |
| `results/` | *not in the repository* – outputs (`.h5ad`, figures) |

Data are not shared in this repository. To use the manual, export cell measurements from QuPath
(*Measure → Export measurements → Cells*) into `Cell_dectection_results/`.
The pipeline was derived from two local notebooks that aren't included here: `Skin_1_analysis.ipynb`
(original pipeline) and `0921_scRNA_IMC.ipynb` (IMC example; source of the area filter and clipping steps).

## Contributing: notebook outputs are stripped on commit

Notebook outputs contain plots of the data, so they are removed when committing while your local
copy keeps them. `.gitattributes` assigns notebooks to the `strip-notebook-output` filter; each clone
must define it once (after `uv sync`):

```bash
git config filter.strip-notebook-output.clean \
  '.venv/bin/jupyter nbconvert --to notebook --ClearOutputPreprocessor.enabled=True --stdin --stdout --log-level=ERROR'
git config filter.strip-notebook-output.smudge cat
```

Without it, outputs would be committed – clear them manually (*Edit → Clear All Outputs*) before committing.

**Exception:** `tma_exploration.ipynb` is the manual, so it is committed *with* its outputs (figures and
tables) so they can be read on GitHub (`tma_exploration.ipynb -filter` in `.gitattributes`). Re-run it
completely before committing so the outputs match the code.

## Running the full pipeline

Run the pipeline on all TMAs/assays:

```bash
uv run python .claude/skills/analyzing-qupath-cells/scripts/run_pipeline.py \
    "Cell_dectection_results/*_cell_detection_TMA*.txt" --out results/tma_pipeline --group-by Parent
```
