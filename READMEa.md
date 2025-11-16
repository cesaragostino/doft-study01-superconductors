# DOFT Pipeline v6 – CSV Structure and Global Patterns (w = 800, p = 7919)

## 1. Contents of the CSV Files (Pipeline Overview)

For the run with **w = 800** and **p = 7919**, the pipeline produces the following files.

## Directory structure

```text
data/         # Materials catalogue (v6), aggregated CSVs, calibration summaries
notebooks/    # Python modules implementing the pipeline (run_all.py + doftstudy/)
paper/        # Manuscript + LaTeX sources + bibliography

Use the supplied Conda specification (Python 3.11 + numpy/pandas/scipy/scikit-learn/matplotlib/tqdm) to mirror the development setup:

```bash
conda env create -f environment.yml
conda activate doft-study01
```
`notebook/run_all.py` delegates to `notebook/doftstudy/run_master_pipeline.py`, which exposes a CLI suitable for both canonical and ad-hoc runs. Example (reproduces the published grid):

```bash
python -m notebook.run_all \
  --input_file data/materials_clusters_real_v6.csv \
  --output_root data/processed \
  --run_all
```

### 1.1 `materials_clusters_real_v6.csv`
252 rows (material + sub_network).
Main fields: `name`, `sub_network`, `category`, `Tc_K`, `Gap_meV`, `ThetaD_K`, `EF_eV`,
`lock_family`, `sub_order`, `notes`.
Dominant sub-network: `single` (207/252).
Remaining rows cover:
* multi-band splits (`band1`, `band2`, `gap1`, `gap2`),
* electron / hole sub-bands (`electron-band`, `hole-band`),
* sigma, pi,
* pressure modes (`0-bar`, `1-bar`, `10-bar`),
* special high-pressure modes (`H1-optic-sigma`, `H2-optic-pi`, `La-acous-sigma`).

### 1.2 `results_calib_w800_p7919.csv`
621 calibrated jumps, each a combination of (`name`, `sub_network`, `jump_desc`).
Key fields:
* `jump_desc` (e.g. `thermal@Tc→gap`, `gap→Debye`, `Debye→E_F`, `sigma-vs-pi`, band vs band, etc.)
* `R_obs`
* `chosen_lock`, `prime_value`
* `err_before`
* `R_corr_eta`, `err_after_eta`
* `R_corr_kappa`, `err_after_kappa`
In this run, `err_after_kappa` is `NaN` everywhere: this is an `η`-only calibration.
`κ` is only used later, at the cluster and fingerprint stages.

### 1.3 `results_cluster_nokappa_w800_p7919.csv` and `results_cluster_kappa_w800_p7919.csv`
175 rows each: the cluster-core subset of jumps.
Same key (`name`, `sub_network`, `jump_desc`, ...) in both files; all 175 entries match one-to-one.
Their `κ`–impact is summarized in `kappa_impact_v6.csv`:
* 169/175 jumps have `Delta_error = 0` within numerical precision.
* Only **6 jumps** show non-zero `Delta_error`, all corresponding to MgB₂ (sigma / pi).
Interpretation: `κ` acts as a local refinement for MgB₂; the rest of the cluster core is effectively invariant under `κ`.

### 1.4 `results_fp_*_full_factorized.csv` (κ / no-κ)
796 rows each (`κ` and `no-κ` versions).
These files carry the full fingerprint information per jump:
* `log_residual_eta`
* `exp_a_2`, `exp_b_3`, `exp_c_5`, `exp_d_7`, `exp_other`
* `p`, `q` (rational lock)
In the `no-κ` file:
* **327** rows have complete factorization in the prime basis (all four exponents defined);
* the remaining rows cannot be expressed cleanly in this simple basis.
In v6, `κ` does **not** change the prime structure:
* For all matched rows, the exponents and `q` are exactly the same in the `κ` and `no-κ` versions.
* `κ` only perturbs the corrected ratios and errors (`R_corr_kappa`, `err_after_kappa`), not the discrete fingerprint itself.

### 1.5 `integer_fingerprint_summary_v6.csv`
19 rows, aggregated by (`category`, `sub_network`) for **integer** locks.
Fields:
* `N` (number of integer fingerprints in that group),
* `exp2_mean`, `exp2_std`,
* `exp3_mean`, `exp3_std`,
* `exp5_mean`, `exp5_std`,
* `exp7_mean`, `exp7_std`.
This is the compact summary used in the integer-fingerprint figures
(per-family prime exponents).

### 1.6 `rational_q_summary_v6.csv` and `kappa_impact_v6.csv`
**`rational_q_summary_v6.csv`:**
* 9 rows, one per superconducting category (Binary, HighPressure, IronBased,
    Molecular, Oxide, HeavyFermion, etc.).
* Fields: `N`, `q_mean`, `q_std`, `q_median`, `q_min`, `q_max`.
* Encodes how each family sits on the rational grid (distribution of denominators `q`).

**`kappa_impact_v6.csv`:**
* 175 rows, one per cluster-core jump.
* Fields: `Delta_error` measuring the change in post-correction error when `κ` is switched on.
* As noted above, only **6 MgB₂** jumps have non-zero `Delta_error`.

---

## Summary: What the Pipeline Encodes
* **Prime-space fingerprint** — exponents on 2, 3, 5, 7.
* **Log-residual layer** — deviations by family via `log_residual_eta`.
* **Rational structure** — `p/q` fractions and their denominators `q`.
* **Calibration layer** — global `Γ`, `η` with bootstrap CIs across different (w, p) settings.

---

## 2. Global Patterns (no-κ, w = 800, p = 7919)
All the following refer to the `no-κ` fingerprints from
`results_fp_nokappa_w800_p7919_full_factorized.csv`.

### 2.1 Log-residuals by family
Grouping `log_residual_eta` by (`category`, `sub_network`):

**Groups with the largest negative means (most compressed relative to the locking grid):**
* `SC_HighPressure / La-acous-sigma`
    mean ≈ −0.29, std ≈ 0.61, n = 3
* `SC_TypeII / single`
    mean ≈ −0.14, std ≈ 0.38, n = 54
* `SC_Binary / gap2`
    mean ≈ −0.12, std ≈ 0.31, n = 6
* `SC_HighPressure / H1-optic-sigma`
    mean ≈ −0.09, std ≈ 0.09, n = 3
* `Superfluid / 0-bar`
    mean ≈ −0.068, std ≈ 0.033, n = 3
* `Superfluid / 10-bar`
    mean ≈ −0.031, std ≈ 0.038, n = 6
* `SC_TypeI / single`
    mean ≈ −0.036, std ≈ 0.068, n = 27
* `SC_Binary / single`
    mean ≈ −0.053, std ≈ 0.16, n = 174
→ These groups tend to fall *below* the rational-locking grid: the observed `R`
is systematically smaller than the ideal prime-locked value.

**Very tight groups (small std, n ≥ 5, well centered near zero):**
* `SC_Binary / pi`
    mean ≈ −0.0024, std ≈ 0.0099, n = 9
* `SC_Molecular / single`
    mean ≈ −0.0030, std ≈ 0.0118, n = 30
* `SC_IronBased / single`
    mean ≈ −0.0023, std ≈ 0.0205, n = 36
* `SC_Oxide / single`
    mean ≈ −0.0051, std ≈ 0.0152, n = 42
* `SC_HighPressure / single`
    mean ≈ −0.0031, std ≈ 0.0382, n = 240
→ These form the **canonical network** where the log-residuals are tightly aligned
with the universal correction law.

### 2.2 Prime structure (no κ)
All statistics below are computed from
`results_fp_nokappa_w800_p7919_full_factorized.csv`.

**Rational denominators `q`**
`q` is defined in **469** rows.
* Range: `q` ∈ {1, 2, 3, 4, 5, 6, 7, 8}.
* Counts:
    * q = 5 → 110
    * q = 8 → 100
    * q = 7 → 74
    * q = 3 → 67
    * q = 4 → 42
    * q = 1 → 40
    * q = 2 → 30
    * q = 6 → 6
* Summary statistics:
    * mean q ≈ **5.06**
    * median q = **5**
→ The typical rational lock lives at **small denominators (3–8)**, with a clear bias
towards **5, 8, 7**. This is consistent with a clean DOFT-like signature: the
natural scale network selects fractions with small denominators.

**Average exponents (327 rows with complete factorization)**
Considering only rows with all four exponents defined:
* ⟨exp₂⟩ ≈ **1.50**
* ⟨exp₃⟩ ≈ **0.65**
* ⟨exp₅⟩ ≈ **0.45**
* ⟨exp₇⟩ ≈ **0.26**

Medians:
* med(exp₂) = **1**
* med(exp₃) = med(exp₅) = med(exp₇) = **0**
→ In blunt terms: most locks are effectively **powers of 2** with discrete
corrections in 3, 5 and 7. The role of 3/5/7 is more of **fine modulation**
than of coarse structure.

**Typical family-level patterns (`integer_fingerprint_summary_v6.csv`)**
Some examples of the average integer fingerprints per group:
* `SC_Binary / pi`:
    rich mixture of 2,3,5,7 with exponents of order unity.
* `SC_TypeI / single`:
    strong weight on 2 and 3, with 5 and 7 present but smaller.
* `SC_TypeII / single`:
    even stronger weight on 2: very “binary”.
* `SC_IronBased / sigma`:
    dominated by 2 and 7; 3 and 5 are almost switched off.
* `HighPressure / H1-optic-sigma`:
    principal axis in 2 and 5.
→ Each macroscopic sub-network ends up with a characteristic “chord” in
{2, 3, 5, 7}.

---

## 3. Effect of κ vs no-κ
### 3.1 Cluster level (175 jumps)
`kappa_impact_v6.csv` summarizes the change in error per jump.
* **169/175** jumps: `Delta_error = 0` within numerical precision.
* **6/175** jumps (all in `MgB₂`, sigma / pi, different thermal and `Debye→E_F` jumps):
    * |Δerror| of order 10⁻³–10⁻².
* **Conclusion**: `κ` only produces appreciable adjustments in MgB₂, one of the most
    complicated systems (multi-gap, strongly coupled). The rest of the cluster
    core is essentially invariant under `κ`.

### 3.2 Full factorization (796 rows)
Comparing `results_fp_kappa_*` vs `results_fp_nokappa_*`:
* For all 796 rows, the exponents
    `exp_a_2`, `exp_b_3`, `exp_c_5`, `exp_d_7`
    and the denominator `q` are **identical**.
* `κ` only modifies `R_corr_kappa` and `err_after_kappa`, not the discrete
    prime-space fingerprint.

### 3.3 Aggregated structure
* At the family level (via `integer_fingerprint_summary_v6.csv` and
    `rational_q_summary_v6.csv`), the means and confidence intervals of the
    exponents and `q` are the same with `κ` and without `κ` in this run.
→ From the perspective of the coarse topology in prime space,
the universal `η` correction captures almost all of the structure;
`κ` refines a few local details without altering the global fingerprint.

