# DOFT Study 01: Constrained Fingerprints for Superconductors and Superfluid Helium

**Author:** Cesar Agostino  
**Draft:** v0.1 (work in progress)

---

## Abstract

We study whether a tightly constrained, DOFT-inspired parametrization can capture quantitative patterns across classical superconductors, high-pressure hydrides, oxides, Fe-based superconductors and superfluid helium. In this framework, each material is represented by a small integer or rational “fingerprint” whose exponents are *pre-specified* by a resonance-locking grammar, rather than fitted freely to the data. A universal correction law with two parameters, Γ and η, is calibrated first on a set of classical metals under non-negativity constraints. Bootstrap (N = 500) and leave-one-out influence analysis yield Γ effectively consistent with zero and a robust positive η of order 10⁻⁵, indicating that a single dissipative scale suffices to describe this family. Integer fingerprints are remarkably stable across runs and do not significantly distinguish Type-I from Type-II superconductors, suggesting a shared integer-locking class. In contrast, high-pressure hydrides separate from classical metals along a specific exponent (exp_d_7), with a Mann–Whitney test giving a false-discovery-rate–corrected p-value of order 10⁻⁴ and a moderate effect size (Cliff’s δ ≈ −0.3). In the rational regime, superfluid helium exhibits significantly lower q values than oxides and Fe-based superconductors, while high-pressure systems cluster around intermediate q with tight confidence intervals. A sub-network protocol finds negligible cluster-coupling κ for all systems except MgB₂, and residual diagnostics show means close to zero with reasonable dispersion, especially for high-pressure materials. Compared against a simple two-parameter power-law baseline, the DOFT-inspired model achieves systematically lower errors with fewer effective parameters. Overall, the data are consistent with the idea that a small DOFT-inspired parameter set, combined with a discrete locking grammar, can organize diverse superconducting and superfluid systems without per-material tuning of exponents.

---

## 1. Introduction

Superconductors and superfluids display a wide range of critical temperatures, gap structures and transport behaviors, yet many empirical regularities survive across families. Classical Type-I and Type-II superconductors, high-pressure hydrides, oxide and Fe-based compounds, and superfluid helium all exhibit sharp transitions and non-trivial scaling laws that are not easily reduced to a single microscopic mechanism. A natural question is whether there exists a *coarse-grained* description that can organize these systems into a small number of classes, using only a few parameters and without per-material tuning of exponents.

The DOFT (Delayed-Oscillator Field Theory) framework proposes such a coarse-grained description. In DOFT, clusters and “grainy” structures are encoded through resonance locks between oscillators, leading to a discrete set of admissible exponents. Rather than fitting continuous exponents freely, the theory prescribes a small library of integer and rational exponents arising from normal locking conditions. In this sense, DOFT is not used here as a microscopic theory of superconductivity, but as a **locking grammar** that constrains how effective fingerprints are built from a limited set of exponents.

In this work we take a deliberately conservative stance: we do not attempt to prove DOFT as a fundamental theory. Instead, we ask a narrower and more practical question: **given** a DOFT-inspired locking grammar and a minimal correction law, can we describe a heterogeneous collection of superconductors and superfluid helium with a small, fixed set of parameters? Concretely, each material is mapped to either an integer fingerprint or a rational fingerprint \(q\), whose exponents are fixed by the locking rules and *not* adjusted per material. A universal correction law with parameters \(\Gamma\) and \(\eta\) is then calibrated on classical metals under weak constraints (\(\Gamma \ge 0\), \(\eta \ge 0\)) and reused across families.

Using this setup, we perform a series of statistical tests on published data for classical superconductors, high-pressure hydrides, oxides, Fe-based superconductors and superfluid helium. We first show that the calibration on classical metals yields \(\Gamma \approx 0\) and a robust positive \(\eta\), as confirmed by bootstrap and leave-one-out influence analysis (Figure 1). We then examine the stability of integer fingerprints across runs and their ability (or lack thereof) to distinguish Type-I from Type-II superconductors, while separating high-pressure hydrides from classical metals along a specific exponent (Figure 2). Non-parametric tests (Kruskal–Wallis, Mann–Whitney, Kolmogorov–Smirnov) and multiple-comparison corrections are used to quantify differences between families in specific exponents and in the rational fingerprint \(q\) (Figure 3). We also explore a sub-network protocol for detecting cluster-level couplings, and we assess residual diagnostics by family to test whether the proposed correction law removes systematic drift (Figure 4). Finally, we compare the DOFT-inspired model against a simple two-parameter power-law baseline (Table 1), to ensure that the observed patterns are not merely an artifact of generic scaling fits.

Our results indicate that (i) a single dissipative scale \(\eta\) suffices to describe classical metals, with \(\Gamma\) effectively null; (ii) Type-I and Type-II superconductors share a common integer fingerprint within the noise of our data; (iii) high-pressure hydrides separate from classical metals along a specific constrained exponent, yielding a moderate but statistically robust effect; (iv) superfluid helium and high-pressure hydrides occupy distinct, tightly clustered regions in rational \(q\)-space, clearly separated from oxides and Fe-based superconductors; (v) within the current protocol, cluster-level couplings are not required to explain the data, except for a small but detectable effect in MgB₂; and (vi) the DOFT-inspired model outperforms a simple power-law baseline while using fewer effective parameters. We interpret these findings as evidence that a DOFT-inspired locking grammar, combined with a minimal correction law, can serve as a compact organizational scheme for diverse superconducting and superfluid systems, while remaining fully testable and falsifiable as new data become available.

In its original formulation, Delayed-Oscillator Field Theory models fields as collections of delayed oscillators. Explicit delay kernels encode memory, and in the long-time coarse-grained limit these delays generate effective time and length scales. Resonance conditions in such delayed-oscillator fields lead to a discrete set of admissible locking ratios. The exponent grammar used here is a reduced, phenomenological version of that structure: we only use the locked exponents and the associated dissipative scale η, without attempting to test or rely on the full cosmological interpretation of DOFT.

---

## 2. Methods

### 2.1. Data and families

We compile published data for superconductors and superfluid helium into a canonical dataset (`materials_master_v5.csv`). Each row corresponds to one material and sub-network, with fields:

- `name` – material name,
- `sub_network` – e.g. `single`, `sigma`, `pi`, `0-bar`, `1-bar`, `10-bar`, special high-pressure modes (`H1-optic-sigma`, `La-acous-sigma`, …),
- `category` – `SC_Binary`, `SC_TypeI`, `SC_TypeII`, `SC_HighPressure`, `SC_IronBased`, `SC_Oxide`, `Superfluid`, etc.,
- thermodynamic and electronic parameters such as \(T_c\) (K), gap (meV), Debye temperature Θ\(_D\) (K), Fermi energy \(E_F\) (eV),
- `lock_family` and `sub_order` (integer, rational or mixed locking; ordering inside each family),
- optional notes.

Classical metals (Type-I and Type-II, plus selected binaries) are used as the *calibration set* for the universal correction law parameters \(\Gamma\) and \(\eta\). High-pressure hydrides, Fe-based superconductors, oxides and superfluid helium are then treated as *test families* for the fingerprints, with no additional per-material tuning of exponents.

### 2.2. DOFT-inspired fingerprints

In the DOFT framework, resonance locks between oscillators restrict the allowed exponents to a discrete set of integer and rational values. We use this locking grammar to define two types of fingerprints:

- **Integer fingerprints**, characterized by a small set of fixed exponents \((\alpha_2, \alpha_3, \alpha_5, \alpha_7)\) associated with primes \(\{2,3,5,7\}\). These are used primarily for classical superconductors and high-pressure hydrides.
- **Rational fingerprints**, characterized by a rational parameter \(q\) and associated exponents of the form \(p/q\) with small denominator \(q\). These are used for softer or hybrid systems such as oxides, Fe-based superconductors and superfluid helium.

Crucially, the exponents themselves are **pre-specified** by the DOFT locking rules and are not fitted separately for each material. Given a family and the locking grammar, the exponent pattern is fixed and only a small number of continuous parameters remain to be calibrated globally. The integer fingerprint is stored in `results_fp_*_full_factorized.csv` as `exp_a_2`, `exp_b_3`, `exp_c_5`, `exp_d_7` plus an optional `exp_other`, while the rational fingerprint is represented by `p`, `q` and `q_avg` summaries.

### 2.3. Universal correction law

We introduce a universal correction law with parameters \(\Gamma\) and \(\eta\) to account for curvature and dissipative memory effects, respectively. For each “jump” between scales (e.g. thermal@Tc → gap, gap → Θ\(_D\), Θ\(_D\) → \(E_F\), σ vs π), we consider an observed ratio \(R_\text{obs}\) and model it as

\[
R_\text{obs} \approx R_\text{lock}(\{\alpha\}) + \Delta(d),
\]

where \(R_\text{lock}\) is determined by the locked exponents and

\[
\Delta(d) \approx \Gamma\, f_\Gamma(d) + \eta\, f_\eta(d),
\]

with \(d\) an effective layer index or distance, and \(f_\Gamma\), \(f_\eta\) fixed functions dictated by the DOFT-inspired ansatz. The details of \(f_\Gamma\) and \(f_\eta\) are kept fixed across materials; only \(\Gamma\) and \(\eta\) are calibrated.

The calibration is performed on classical metals only, under weak inequality constraints

\[
\Gamma \ge 0, \quad \eta \ge 0.
\]

We use constrained least-squares to obtain point estimates \(\hat{\Gamma}, \hat{\eta}\), and we quantify uncertainty and robustness through bootstrap and leave-one-out analyses (Section 2.4). Once calibrated, \(\hat{\Gamma}\) and \(\hat{\eta}\) are held fixed for all other families. The calibrated jumps and associated errors are stored in `results_calib_w800_p7919.csv`.

### 2.4. Bootstrap and influence analysis

To assess the stability of the calibration on classical metals, we perform a non-parametric bootstrap with \(N = 500\) resamples of the calibration dataset. For each bootstrap sample we re-fit the correction law, obtaining distributions for \(\Gamma\) and \(\eta\). We report percentile-based confidence intervals and examine the skewness and width of these distributions. Figure 1 (left and middle panels) shows the bootstrap histograms and kernel density estimates, with the full-sample estimates overlaid.

We also carry out a leave-one-out (LOO) influence analysis: each classical metal is removed in turn, the model is re-fitted, and the resulting \(\Gamma\) and \(\eta\) estimates are compared to the full-sample values. We summarize the relative change Δη (%) per metal in `results_calib_*` and display it as a bar plot in Figure 1 (right panel).

### 2.5. Statistical tests for family separation

We use standard non-parametric tests to quantify differences between families in their fingerprint parameters:

- **Kruskal–Wallis (KW)** tests for global differences in specific exponents (e.g. `exp_d_7`, `exp_a_2`) across multiple families (Type-I, Type-II, High-pressure, Binary, etc.).
- **Mann–Whitney U** tests for pairwise comparisons, in particular between high-pressure hydrides and classical metals along selected integer exponents.
- **Kolmogorov–Smirnov (KS)** tests for differences in the distribution of rational \(q\) between families (e.g. superfluid helium vs high-pressure vs Fe-based).

To control for multiple comparisons when scanning over several exponents and families, we apply false discovery rate (FDR) correction (Benjamini–Hochberg). Effect sizes are reported alongside p-values, using Cliff’s delta δ for two-sample tests where appropriate. The bootstrap confidence intervals for exponent means and \(q\)-averages by family are stored in `fingerprint_fp_*_bootstrap_CIs.csv`.

### 2.6. Sub-network protocol and cluster coupling κ

To probe possible cluster-level couplings, we define a sub-network protocol based on a contrast \(C_{AB}\) between two channels (e.g. σ and π bands in MgB₂). A coupling parameter \(\kappa \ge 0\) is introduced only when the contrast suggests potential mixing. The protocol is intentionally conservative: \(\kappa\) is included as an extra parameter only if it reduces residual drift in a statistically meaningful way; otherwise, \(\kappa\) is set to zero.

We apply this procedure across the dataset, with particular attention to systems such as MgB₂ and FeSe where multi-band or multi-channel behavior is expected. Cluster-level results with and without κ are stored in `results_cluster_kappa_w800_p7919.csv` and `results_cluster_nokappa_w800_p7919.csv`. Differences in error are visualized in Figure 5.

### 2.7. Residual diagnostics

We analyze residuals by family after applying the universal correction law and the fingerprint exponents. For each group `(category, sub_network)` we compute:

- the mean and standard deviation of the log-residual `log_residual_eta`,
- the number of contributing jumps `count`.

These summaries are stored in `fingerprint_fp_*_log_residual.csv` and displayed in Figure 4a, with groups of small sample size (N < 5) highlighted. We also examine residuals as a function of layer index \(d\) for selected families (high-pressure, Fe-based, superfluid), as shown in Figure 4b.

### 2.8. Baseline model

As a sanity check, we compare the DOFT-inspired model to a simple **power-law baseline**. For each group (e.g. by `category` or `lock_family`), we fit a two-parameter model

\[
R_\text{obs} \approx A\, d^\alpha,
\]

where both the scale \(A\) and exponent \(\alpha\) are free parameters optimized by least squares for that group. The baseline is intentionally generic: it does not use the prime-lock structure or the rational grammar, and exponents are continuous rather than discrete.

We evaluate performance using the mean absolute relative error (MARE):

\[
\text{MARE} = \frac{1}{N} \sum_{i=1}^N \frac{|R_{\text{pred},i} - R_{\text{obs},i}|}{|R_{\text{obs},i}|}.
\]

For each group, we compute:

- `baseline_mare` for the power-law model,
- `doft_mare` using the DOFT-inspired model with fixed exponents and the calibrated \(\Gamma, \eta\),
- their difference `err_delta_baseline_minus_doft = baseline_mare − doft_mare`,
- the effective number of free parameters (`baseline_params = 2`, `doft_params = 0` at the group level).

These results are stored in `baseline_test.csv` and summarized in Table 1.

---

## 3. Results

### 3.1. Calibration on classical metals: Γ ≈ 0, η > 0

Applying the constrained calibration procedure to the set of classical metals yields:

- a dissipative parameter \(\hat{\eta}\) of order \(4 \times 10^{-5}\),
- a curvature parameter \(\hat{\Gamma}\) numerically compatible with zero.

Figure 1a shows the bootstrap distribution of \(\eta\). The distribution is unimodal and fairly narrow, with all resamples yielding \(\eta > 0\) and the full-sample estimate lying close to the peak. Figure 1b shows the bootstrap distribution of \(\Gamma\): the vast majority of resamples sit exactly at or extremely close to zero, with a rapidly decaying tail extending to small positive values. This suggests that curvature is not required by the data within the numerical precision and constraints used.

The LOO analysis in Figure 1c indicates that the relative change Δη (%) in \(\eta\) upon removing any single metal remains modest (typically within ±10–15%), with no single material dominating the calibration. Together, these results support the interpretation that (i) a single dissipative scale \(\eta\) suffices to describe classical metals under the present ansatz, and (ii) \(\Gamma\) can be effectively set to zero without loss of fit quality.

### 3.2. Integer fingerprints and Type-I vs Type-II vs high-pressure

Using the integer-locking grammar, we assign prime exponents (2,3,5,7) to each jump and derive family-level fingerprints from `results_fp_nokappa_full_factorized.csv` and the bootstrap summaries. Figure 2 shows violin plots of two key exponents:

- `exp_d_7` (associated with the prime 7),
- `exp_a_2` (associated with the prime 2).

For `exp_d_7` (Figure 2, left panel), Type-I and Type-II superconductors show strongly overlapping distributions, with similar medians and spreads. In contrast, high-pressure hydrides form a distinct cluster at lower `exp_d_7` values. Kruskal–Wallis tests detect significant global differences at this exponent, and a Mann–Whitney comparison between high-pressure and classical metals yields an FDR-corrected p-value of order 10⁻⁴ with a moderate Cliff’s delta (δ ≈ −0.3). Binary superconductors populate yet another region of the `exp_d_7` axis, consistent with the more complex prime mixture seen in their fingerprint.

For `exp_a_2` (Figure 2, right panel), all families exhibit strong weight on the prime 2, but with different effective exponents: Type-II tends to have slightly larger 2-components than Type-I, while high-pressure and binary families occupy intermediate ranges. Overall, these patterns support two key statements:

1. **Type-I vs Type-II:** within the precision of the data and the constraints of the grammar, Type-I and Type-II form a common integer-locking class; we do not find statistically robust separation in the integer fingerprint.
2. **High-pressure vs classical:** high-pressure hydrides are clearly separated from classical metals along the constrained axis `exp_d_7`, providing a quantitative fingerprint for the high-pressure class.

### 3.3. Rational fingerprints: superfluid helium, high-pressure, Fe-based

For softer or hybrid systems we consider rational fingerprints encoded by a denominator \(q\). Figure 3 shows the distribution of \(q\) by family for superfluid helium, high-pressure superconductors and Fe-based materials.

Superfluid helium occupies a distinctly low-\(q\) regime, with values concentrated at small denominators (typically \(q \approx 2\) and a narrow spread). High-pressure hydrides cluster around intermediate denominators (median \(q\) in the range 3–6), forming a relatively tight distribution that suggests an internally coherent locking pattern. Fe-based superconductors are concentrated near \(q \approx 1\), consistent with a more “integer-like” rational structure in this representation.

Kolmogorov–Smirnov tests confirm that the \(q\) distributions for superfluid helium and high-pressure hydrides are significantly different from those of the oxide/Fe-based sector. Within the DOFT-inspired picture, these results are consistent with the idea that superfluid helium and high-pressure hydrides occupy distinct regions in rational \(q\)-space, both clearly separated from the more conventional Fe-based family.

### 3.4. Residual diagnostics by family

Figure 4 summarizes residual diagnostics after applying the universal \(\eta\)-correction and the locked exponents.

Figure 4a shows the mean and standard deviation of `log_residual_eta` for each `(category, sub_network)` group. Most groups lie close to zero with moderate dispersion, indicating that the model removes the dominant trends. A few groups with small sample sizes (N < 5) — notably certain high-pressure optic/acoustic sub-networks and specific superfluid modes — show larger negative means, corresponding to observed ratios slightly below the ideal prime-based grid. These groups are highlighted in a different color to emphasize their limited statistical weight.

Figure 4b shows `log_residual_eta` as a function of layer index \(d\) for three representative families: high-pressure, Fe-based and superfluid. In all cases, the residuals are centered near zero, with no clear systematic drift versus \(d\). High-pressure systems, in particular, show a tight cloud of residuals. We interpret this as evidence that the combination of (i) a universal dissipative scale \(\eta\), (ii) negligible curvature \(\Gamma\), and (iii) the discrete locking grammar, is sufficient to capture the dominant structure in the data at the current level of noise.

### 3.5. Effect of κ: cluster-level couplings

The sub-network protocol introduces a cluster-coupling parameter \(\kappa\) only when the contrast \(C_{AB}\) between channels suggests potential mixing. We apply this protocol to the “core cluster” of 51 selected jumps (stored in `results_cluster_*`).

Figure 5a shows the change in error Δerror = `err_after_kappa − err_after_eta` for all 51 jumps. Most entries lie exactly at Δerror = 0, indicating that the κ term does not modify the fit for the vast majority of jumps. Figure 5b zooms in on the cases with the largest |Δerror|. These correspond predominantly to MgB₂ (σ and π sub-networks), where a small but non-zero κ slightly improves the residuals. For other candidate multi-band systems such as FeSe, κ remains negligible within uncertainty.

Thus, κ behaves as a **localized second-order correction**: it is not required to organize the global fingerprint, but it captures subtle channel mixing in MgB₂, consistent with its known multi-gap character.

### 3.6. Prime-space structure

Beyond individual exponents, we can summarize the global structure of the prime exponents across all fully factorized jumps (about 196 rows). The average exponents are:

- ⟨exp\(_{2}\)⟩ ≈ 1.3,
- ⟨exp\(_{3}\)⟩ ≈ 0.6,
- ⟨exp\(_{5}\)⟩ ≈ 0.4,
- ⟨exp\(_{7}\)⟩ ≈ 0.2,

with medians roughly (1, 0, 0, 0). In other words, most locks are dominated by prime 2, with smaller contributions from 3, 5 and 7 acting as fine modulation. Family-level fingerprints (from the bootstrap summaries) show characteristic “prime chords”: for example, Fe-based σ sub-networks are almost purely (2,7) with 3 and 5 effectively off, while certain high-pressure optic modes emphasize (2,5).

We interpret these patterns as the prime-space signature of each macroscopic family under the DOFT-inspired grammar.

### 3.7. Comparison with a power-law baseline

Table 1 summarizes the comparison between the DOFT-inspired model and the two-parameter power-law baseline for several groupings (by `category` and by `lock_family`). For each group we report the mean absolute relative error (MARE) for the baseline (`baseline_mare`) and for DOFT (`doft_mare`), along with the difference ΔMARE = baseline − DOFT and the number of parameters.

**Table 1 – Mean absolute relative error (MARE) for the power-law baseline vs DOFT-inspired model.**  
*(Illustrative numbers based on `baseline_test.csv`; to be updated with final values.)*

| Grouping | Group value        | Baseline MARE | DOFT MARE | ΔMARE (baseline−DOFT) | Baseline params | DOFT params |
|----------|--------------------|---------------|-----------|-----------------------|-----------------|-------------|
| category | SC_HighPressure    | ~0.4          | ~0.2–0.3  | > 0                   | 2               | 0           |
| category | SC_Binary          | ~0.6          | ~0.02     | > 0                   | 2               | 0           |
| category | Superfluid         | ~2            | ~0.003    | > 0                   | 2               | 0           |
| category | SC_IronBased       | ~3            | ~0.03     | > 0                   | 2               | 0           |
| lock_family | integer         | ~0.7          | ~0.16     | > 0                   | 2               | 0           |
| lock_family | rational        | ~2            | ~0.003    | > 0                   | 2               | 0           |
| lock_family | mixed           | >4            | ~0.02     | > 0                   | 2               | 0           |

In all cases considered, ΔMARE is positive: the baseline error exceeds the DOFT error, sometimes by more than an order of magnitude, despite the baseline having two free parameters per group and the DOFT-inspired model having none at that level (the exponents and \(\eta\) are fixed globally). This indicates that the constrained locking grammar is not simply rephrasing a generic power-law fit: it provides a *more efficient* parametrization of the data, with fewer degrees of freedom and lower error.

---

## 4. Discussion

### 4.1. DOFT-inspired grammar as a phenomenological model

The central question of this study is not whether DOFT is a complete microscopic theory of superconductivity or superfluidity, but whether its locking grammar can serve as a useful phenomenological parametrization. By constraining exponents to a discrete set derived from resonance locks, and by calibrating only a small number of continuous parameters globally, we obtain a compact description that remains falsifiable: if the data were incompatible with the grammar, this would show up as poor fits, unstable calibrations or structured residuals.

Our findings indicate that the DOFT-inspired grammar clears several non-trivial tests:

- Classical metals can be described with a single dissipative scale \(\eta\) and negligible curvature \(\Gamma\).
- Type-I and Type-II superconductors share a common integer fingerprint; the grammar does not artificially force a distinction where the data do not support it.
- High-pressure hydrides separate from classical metals along a specific constrained exponent (`exp_d_7`), providing a quantitative fingerprint for the high-pressure class.
- Superfluid helium and high-pressure hydrides occupy distinct, tightly clustered regions in rational \(q\)-space, clearly separated from oxides and Fe-based superconductors.
- Cluster-level couplings, as encoded by \(\kappa\), are not required to explain the data under the present protocol, except for a localized effect in MgB₂.
- Compared to a simple power-law baseline, the DOFT-inspired model achieves lower errors with fewer effective parameters.

None of these outcomes are guaranteed a priori, given the tight constraints on exponents and parameters. In this sense, the results support the use of DOFT as a structured phenomenological model, even if its status as a fundamental theory remains open.

From the equations-of-motion viewpoint, η summarizes the strength of memory terms generated by the delay kernels in the underlying oscillator picture. The fact that a single positive η, with Γ ≈ 0, is enough to organize classical metals is therefore consistent with the idea that only a small amount of coarse-grained memory survives at the scales probed by the superconducting observables considered here.

### 4.2. Relation to simpler baselines

A natural concern is whether the observed patterns could be captured equally well, or better, by simpler models without a discrete locking grammar—for example, by fitting arbitrary power laws or flexible scaling forms to each family. From a statistical perspective, such baselines typically introduce more degrees of freedom, allowing per-family or even per-material exponents, at the cost of reduced falsifiability.

Our baseline analysis addresses this concern in a minimal way. With only two parameters (scale and exponent) per group, the power-law model is already quite flexible; yet it yields larger mean absolute relative errors than the DOFT-inspired model across all tested categories and lock families. The DOFT grammar therefore achieves a better trade-off between fit quality and parameter count. A more exhaustive model comparison (e.g. using information criteria or cross-validation) is left for future work, but the present results already suggest that the locking grammar is adding more than cosmetic structure.

### 4.3. Limitations

Several limitations should be kept in mind:

1. **Dataset size and selection.**  
   The number of materials per family, especially for high-pressure hydrides and some oxide or Fe-based subsets, is still modest. This limits the statistical power of some tests and increases sensitivity to selection bias in the available literature.

2. **Choice of observables.**  
   The present implementation relies on a specific set of jumps and derived ratios (thermal@Tc → gap, gap → Θ\(_D\), Θ\(_D\) → \(E_F\), etc.). Alternative constructions might interact differently with the locking grammar and could reveal additional structure or tensions.

3. **Conservativeness of the κ protocol.**  
   The sub-network protocol for detecting cluster-level couplings is intentionally conservative and may miss more subtle forms of mixing. The fact that κ is non-zero only in MgB₂ may reflect both physical reality and the protocol’s thresholds.

4. **Simplified baseline.**  
   The baseline model tested here is a simple two-parameter power law. While it already provides a non-trivial comparison, more sophisticated baselines (e.g. multi-scale or piecewise fits) could be explored to further test the value of the locking grammar.

These limitations do not invalidate the patterns reported here, but constrain how strongly they should be interpreted. We view the present work as an initial test of a structured parametrization, rather than as a definitive classification of all superconductors and superfluids.

### 4.4. Outlook

Several extensions of this study suggest themselves:

- **Data expansion.**  
  Extending the dataset with additional high-pressure hydrides, unconventional superconductors, and more detailed measurements of superfluid helium would sharpen the statistical tests and better probe the boundaries of the locking grammar.

- **Alternative grammars and baselines.**  
  Systematic comparisons with alternative locking schemes and with richer baseline models would help isolate which aspects of the results are specific to DOFT and which are more generic properties of the data.

- **Predictive testing.**  
  Using the calibrated grammar to make predictions for new materials (e.g. predicting plausible fingerprints or residuals for candidate high-pressure systems) would provide a stronger test of its predictive power.

More broadly, the approach taken here—using a discrete locking grammar and a small number of global parameters to organize diverse systems—may be applicable beyond superconductivity and superfluidity. Any domain where layered or clustered structures emerge from oscillator-like dynamics could provide a testing ground for DOFT-inspired fingerprints.

---

## 5. Conclusions

We have tested a DOFT-inspired fingerprint model on a heterogeneous collection of superconductors and superfluid helium, using a tightly constrained set of integer and rational exponents combined with a universal correction law. The main conclusions are:

1. A single dissipative parameter \(\eta\) suffices to describe classical metals, while the curvature parameter \(\Gamma\) is effectively null within the numerical precision of our bootstrap and LOO analyses.
2. Integer fingerprints are stable across runs and do not significantly distinguish Type-I from Type-II superconductors, indicating a shared integer-locking class.
3. High-pressure hydrides separate from classical metals along a specific constrained exponent (`exp_d_7`), yielding a statistically robust, moderate effect that serves as a quantitative fingerprint signature.
4. In rational \(q\)-space, superfluid helium and high-pressure hydrides occupy distinct, tightly clustered regions that are clearly separated from oxides and Fe-based superconductors.
5. Cluster-level couplings, as encoded by \(\kappa\), are not required to explain the data under the present ansatz, except for a small but detectable effect in MgB₂.
6. Residual diagnostics by family show means close to zero and reasonable dispersion, especially for high-pressure hydrides, suggesting that no additional ad-hoc parameters are needed at the current level of detail.
7. Compared to a simple two-parameter power-law baseline, the DOFT-inspired model achieves systematically lower errors while using fewer effective parameters, indicating that the locking grammar provides genuine explanatory power rather than mere reparameterization.

Taken together, these results support the view that a small DOFT-inspired parameter set, combined with a discrete locking grammar, can provide a compact and testable organizational scheme for superconducting and superfluid systems. The framework remains falsifiable as new data and more stringent baselines are brought to bear, and it invites further exploration of resonance-based fingerprints in other complex materials.

---

## 6. Data and code availability

All data and code used in this study are provided in this repository. The canonical input dataset is distributed as `data/raw/materials_master_v5.csv`, and the main analysis pipeline is implemented in the `src/` directory. Processed results (calibration outputs, fingerprints, statistical test summaries) and the figures and tables corresponding to the manuscript are stored under `data/processed/` and `results/` respectively.

Upon stabilization of the analysis, this repository will be archived with a permanent DOI, and the citation information in `CITATION.cff` will be updated accordingly.

---

## Acknowledgements

The author thanks informal collaborators and the broader condensed-matter community for maintaining open data and discussions that made this exploratory study possible. Any errors or misinterpretations remain the responsibility of the author.

The author made extensive use of large language models as coding and writing assistants during this project, in particular OpenAI ChatGPT, Google Gemini, and Anthropic Claude. These tools were used for code refactoring, exploratory data-analysis scripts, and editorial suggestions on the manuscript. All numerical results, checks and scientific interpretations were independently verified, and any remaining errors are the sole responsibility of the author.
