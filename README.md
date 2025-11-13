# DOFT Study 01 – Superconductors and Superfluid Helium

This repository contains the code and data for DOFT Study 01, where we test a constrained,
DOFT-inspired fingerprint model on classical superconductors, high-pressure hydrides,
oxides, Fe-based superconductors and superfluid helium.

doft-study01-superconductors/
├─ README.md
├─ LICENSE
├─ CITATION.cff          # opcional pero muy bueno
├─ environment.yml       # o requirements.txt
├─ data/
│  ├─ raw/
│  │  └─ materials_master_v5.csv
│  └─ processed/
│     ├─ results_calib_*.csv
│     ├─ fingerprint_*.csv
│     └─ results_cluster_*.csv
├─ src/
│  ├─ doftstudy/
│  │  ├─ __init__.py
│  │  ├─ calibration.py
│  │  ├─ fingerprints.py
│  │  ├─ stats_tests.py
│  │  └─ plotting.py
│  └─ run_all.py
├─ notebooks/
│  ├─ 01_exploration.ipynb
│  └─ 02_checks.ipynb
├─ results/
│  ├─ figures/
│  │  ├─ fig1_calibration_eta_gamma.png
│  │  ├─ fig2_integer_fingerprints.png
│  │  └─ ...
│  └─ tables/
│     ├─ table1_calibration_summary.csv
│     ├─ table2_stats_tests.csv
│     └─ ...
└─ paper/
   ├─ draft.tex (o .md)
   └─ biblio.bib
