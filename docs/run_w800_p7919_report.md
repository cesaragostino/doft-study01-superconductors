# DOFT Pipeline – CSV Structure and Global Patterns (w = 800, p = 7919)

## 1. Contents of the CSV Files (Pipeline Overview)

For the run with **w = 800**, **p = 7919**, the pipeline produces the following files.

---

### 1.1 `materials_clusters_real_v5.csv`

* **140 rows** (material + sub_network).
* Main fields: `name`, `sub_network`, `category`, `Tc_K`, `Gap_meV`, `ThetaD_K`, `EF_eV`, `lock_family`, `sub_order`, `notes`.
* Dominant sub-network: **single** (127/140).
* Remaining rows: `sigma`, `pi`, `0/1/10-bar`, and high-pressure special modes (`H1-optic-sigma`, etc.).

---

### 1.2 `results_calib_w800_p7919.csv`

* **381 calibrated jumps**, each a combination of (material, sub_network, jump_desc).
* Key fields:

  * `jump_desc` (e.g., thermal@Tc→gap, gap→Debye, Debye→E_F, sigma-vs-pi)
  * `R_obs`
  * `chosen_lock`, `prime_value`
  * `err_before`
  * `R_corr_eta`, `err_after_eta`
* `err_after_kappa` = **NaN** → this is an **η-only** calibration.

---

### 1.3 `results_cluster_nokappa_w800_p7919.csv` and `results_cluster_kappa_w800_p7919.csv`

* **51 rows** each.
* Same key rows in both files; **51 matched** entries.
* Only **6 rows differ**, all corresponding to **MgB₂ (sigma & pi)**.
* κ modifies `R_corr_kappa` and `err_after_kappa` slightly.

---

### 1.4 `results_fp_*_full_factorized.csv` (κ / no-κ)

* **432 rows** each.
* Additional fingerprint info:

  * `log_residual_eta`
  * `exp_a_2`, `exp_b_3`, `exp_c_5`, `exp_d_7`, `exp_other`
  * `p`, `q`
* ~196 rows have complete factorization; the rest cannot be expressed cleanly in the prime basis.

---

### 1.5 `fingerprint_fp_*_log_residual.csv`

* **23 rows** aggregated by (category, sub_network).
* Stores mean, std, and count of `log_residual_eta`.

---

### 1.6 `fingerprint_fp_*_bootstrap_CIs.csv`

* **58 rows**.
* Bootstrap summaries by (family_type, group, lock_type, metric).
* Metrics: `exp_a_2`, `exp_b_3`, `exp_c_5`, `exp_d_7`, `q_avg` with full CI info.

---

## Summary: What the Pipeline Encodes

1. **Prime-space fingerprint** — exponents on 2, 3, 5, 7.
2. **Log-residual layer** — deviations by family.
3. **Rational structure** — p/q fractions and q_avg.

---

# 2. Global Patterns (no-κ, w = 800, p = 7919)

## 2.1 Log-residuals by Family (no κ)

From `fingerprint_fp_nokappa_log_residual.csv`:

Groups most displaced (mean < 0):

* **SC_HighPressure / La-acous-sigma**: mean ≈ −0.29 (n = 3)
* **Superfluid / 0-bar**: mean ≈ −0.1215 (n = 3)
* **SC_HighPressure / H1-optic-sigma**: mean ≈ −0.0942 (n = 3)
* **SC_TypeII / single**: mean ≈ −0.0657 (n = 27)
* **Superfluid / 10-bar**: mean ≈ −0.0545 (n = 6)
* **SC_TypeI / single**: mean ≈ −0.0398 (n = 27)

Interpretation:

* These groups tend to fall **below** the rational-locking grid → observed R is slightly smaller than ideal.
* Superfluids & some HP modes appear **compressed**.

Tight groups (small std, n ≥ 5):

* SC_Binary / pi
* SC_Molecular / single
* SC_IronBased / single
* SC_HighPressure / single

These form the **canonical network** of clean residual alignment.

---

## 2.2 Prime Structure (no κ)

From `results_fp_nokappa_full_factorized.csv`:

### Rational denominators q

* Defined in **236 rows**.
* Distribution peaks: **q = 5, 7, 8**, then 3–4; few at 1,2,6.
* mean q ≈ 5.22, median = 5.

→ The system prefers **small denominators** (3–8). Clean DOFT signature.

### Average exponents (over 196 clean rows)

* ⟨exp₂⟩ ≈ 1.28
* ⟨exp₃⟩ ≈ 0.58
* ⟨exp₅⟩ ≈ 0.39
* ⟨exp₇⟩ ≈ 0.24

Medians:

* exp₂ median = 1
* exp₃, exp₅, exp₇ medians = 0

→ Prime 2 dominates; primes 3,5,7 act as modulation.

### Examples (bootstrap means)

* **SC_Binary / pi**: rich 2,3,7; mild 5.
* **SC_Binary / sigma**: balanced 2,3,5,7.
* **SC_TypeI / single**: strong 2,3; weaker 5,7.
* **SC_TypeII / single**: very strong 2.
* **SC_IronBased / sigma**: almost pure (2,7).
* **SC_HighPressure / H1-optic-sigma**: (2,5)-dominated.

---

# 3. Effect of κ vs no-κ

## 3.1 Cluster Level (51 jumps)

* 45/51 rows identical.
* All 6 differences = **MgB₂**.
* κ only affects borderline cases; shifts are small (~1e−3–1e−2).

## 3.2 Full Factorization (432 rows)

* 426 rows identical.
* 12 rows differ (6 κ only / 6 no-κ only), again **MgB₂**.
* Exponents on 2,3,5,7 are unchanged for matched rows.

## 3.3 Aggregated Log-residuals

* Means, stds, counts identical → κ does **not** change family-level residual structure.

## 3.4 Bootstrap CIs (exponents & q_avg)

Global mean shifts (κ − no-κ):

* exp₂: +0.021
* exp₃: −0.0036
* exp₅: +0.0009
* exp₇: −0.0019
* q_avg: −0.0076 (with local shifts up to ±0.25)

Example (SC_Binary / pi):

* exp₂: 1.243 → 1.310 (↑)
* exp₃: 0.754 → 0.7385 (↓)
* exp₅: 0.2545 → 0.238 (↓)
* exp₇: 0.5085 → 0.4995 (↓)

Interpretation:

* κ applies **second-order corrections**:

  * increases exp₂ slightly
  * slightly perturbs q_avg in some groups
  * preserves the coarse DOFT topology
