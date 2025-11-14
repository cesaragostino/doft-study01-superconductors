# DOFT Study 01 – Superconductors and Superfluid Helium

This repository contains the code and data for **DOFT Study 01**, where we test a constrained,
DOFT-inspired fingerprint model on:

* classical superconductors
* high-pressure hydrides
* oxides
* Fe-based superconductors
* superfluid helium

---

## Top-level Repository Layout

```
doft-study01-superconductors/
  README.md        # Main project documentation
  REFERENCES.md    # Bibliography and external sources
  GOVERNANCE.md    # Project governance, roles, and contribution rules
  CITATION.cff     # Citation metadata for GitHub/Zenodo
  ...              # Additional project files
```

---

## Repository Structure

```
data/
  raw/        # Canonical input data (materials, families, metadata)
  processed/  # Derived CSVs (calibration, fingerprints, clustering)

src/          # Python modules implementing the full pipeline
docs/         # Documentation and run report
notebooks/    # Optional exploration (not required for main results)
results/      # Figures and tables corresponding to the main paper
paper/        # Draft manuscript + bibliography
```

---

## How to Reproduce the Main Results

### 1. Create the environment

Use the supplied Conda specification (Python 3.11 + numpy/pandas/scipy/scikit-learn/matplotlib/tqdm) to mirror the development setup:

```bash
conda env create -f environment.yml
conda activate doft-study01
```

### 2. Run the master pipeline

`src/run_all.py` delegates to `src/doftstudy/run_master_pipeline.py`, which exposes a CLI suitable for both canonical and ad-hoc runs. Example (reproduces the published grid):

```bash
python -m src.run_all \
  --input_file data/raw/materials_clusters_real_v6.csv \
  --output_root data/processed \
  --run_all
```

Common flags:

* `--run_all` – evaluates every predefined combo of winsor caps, prime_max values and jitter settings.
* `--input_file` – override the source CSV (default: `data/raw/materials_clusters_real_v6.csv`).
* `--output_root` – where `results_w<...>_p<...>` folders are created (defaults to the script directory).
* `-w/--winsor_values`, `-p/--prime_values`, `-j/--jitter_values`, `-n/--n_runs` – manual control when `--run_all` is omitted.

Each robustness job produces a structured bundle:

```
results_wXXX_pYYYY/
  calib/         # calibration CSVs + summaries
  cluster/       # cluster diagnostics (κ/no-κ)
  fingerprint/   # DOFT fingerprints, residual logs, bootstrap CIs
  baseline/      # baseline vs DOFT comparison tables
  figures/       # PNG exports (Fig. 1–5)
```

The sensitivity (jitter) stage reuses the `doft_config.json` emitted under the corresponding `calib/` directory.

### 3. Figures and tables used in the paper

Every `results_wXXX_pYYYY/figures/` directory contains the exact PNGs referenced in the manuscript:

| File name                      | Description                             |
| ------------------------------ | --------------------------------------- |
| `fig01_calibration.png`        | Fig. 1 – Calibration (bootstrap + LOO)  |
| `fig02_integer_fingerprint.png`| Fig. 2 – Integer fingerprint exponents  |
| `fig03_rational_q.png`         | Fig. 3 – Rational denominators by family|
| `fig04_residuals.png`          | Fig. 4 – Residual diagnostics           |
| `fig05_kappa_vs_no_kappa.png`  | Fig. 5 – κ vs. no-κ comparison          |

Need to refresh the plots for a specific tag? Run the helper:

```bash
python src/doftstudy/generate_figures.py --base_tag w600_p10000
```

Tables (e.g., `baseline/baseline_summary_<tag>.csv`, calibration summaries, bootstrap stats) live alongside the CSV outputs and can be imported directly into the paper or supplement.

---

## Data

The canonical dataset is:

```
data/raw/materials_clusters_real_v6.csv
```

Each row corresponds to one material.

Key columns include:

* `name` – material name
* `category` – Type-I, Type-II, high-pressure, oxide, Fe-based, superfluid, ...
* `sub_network` – sigma, pi, single, optic/acoustic modes (used for contrasts)
* `Tc_K`, `ThetaD_K`, `EF_eV`, `Gap_meV` – anchor observables for each material
* `lock_family` – preferred locking grammar (`integer`, `rational`, `mixed`)

---

## License & Citation

**License:** MIT

If you use this code or data, please cite:

```bibtex
# TODO: add your paper / preprint reference here
```

---
