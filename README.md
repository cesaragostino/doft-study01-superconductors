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
notebooks/    # Optional exploration (not required for main results)
results/      # Figures and tables corresponding to the main paper
paper/        # Draft manuscript + bibliography
```

---

## How to Reproduce the Main Results

### 1. Create the environment

```bash
conda env create -f environment.yml
conda activate doft-study01
```

### 2. Run the full pipeline

```bash
python -m src.run_all
```

This will:

* read `data/raw/materials_master_v5.csv`
* perform the Γ and η calibration on classical metals
* compute integer & rational DOFT fingerprints
* run all statistical tests (bootstrap, LOO, MWU, KS, FDR, etc.)
* generate:

```
data/processed/     # processed CSVs
results/figures/    # all figures
results/tables/     # summary tables
```

### 3. Figures and tables used in the paper

| File                                             | Description               |
| ------------------------------------------------ | ------------------------- |
| `results/figures/fig1_calibration_eta_gamma.png` | Figure 1                  |
| `results/figures/fig2_integer_families.png`      | Figure 2                  |
| `results/tables/table1_calibration_summary.csv`  | Table 1                   |
| `results/tables/table2_stats_tests.csv`          | Table 2                   |
| `...`                                            | additional figures/tables |

---

## Data

The canonical dataset is:

```
data/raw/materials_master_v5.csv
```

Each row corresponds to one material.

Key columns include:

* `name` – material name
* `family` – classical, HP, oxide, Fe-based, superfluid, ...
* `Tc` – critical temperature (K)
* `...` – include the fields actually used in your analysis

---

## License & Citation

**License:** MIT

If you use this code or data, please cite:

```bibtex
# TODO: add your paper / preprint reference here
```

---
