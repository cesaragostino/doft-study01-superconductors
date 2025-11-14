# DOFT Study 01: Constrained Fingerprints for Superconductors and Superfluid Helium

**Abstract**  
We study whether a tightly constrained, DOFT-inspired parametrization can capture quantitative patterns across classical superconductors, high-pressure hydrides, oxides, Fe-based superconductors and superfluid helium. In this framework, each material is represented by a small integer or rational “fingerprint” whose exponents are *pre-specified* by a resonance-locking grammar, rather than fitted freely to the data. A universal correction law with two parameters, Γ and η, is calibrated first on a set of classical metals under non-negativity constraints. Bootstrap (N = 500) and leave-one-out influence analysis yield Γ ≈ O(10⁻¹⁷) (effectively null curvature) and a robust positive η ≈ 4×10⁻⁵, indicating that a single dissipative scale suffices to describe this family. Integer fingerprints are remarkably stable across runs and do not significantly distinguish Type-I from Type-II superconductors, suggesting a shared “integer-locking” class. In contrast, high-pressure hydrides separate from classical metals along a specific exponent (exp_d_7), with a Mann–Whitney test giving an FDR-corrected p ≈ 10⁻⁴ and a moderate effect size (Cliff’s d ≈ −0.30), which provides a quantitative signature for the high-pressure class. In the rational regime, superfluid helium exhibits significantly lower q values (≈ 2–2.3) than oxides and Fe-based superconductors (Kolmogorov–Smirnov p ≪ 0.01), while high-pressure systems cluster around q ≈ 5.85 with tight confidence intervals. A sub-network protocol finds negligible cluster-coupling κ for all systems except MgB₂, and residual diagnostics show means close to zero with reasonable dispersion, especially for high-pressure materials. Overall, the data are consistent with the idea that a small DOFT-inspired parameter set, combined with a discrete locking grammar, can organize diverse superconducting and superfluid systems without per-material tuning of exponents.

## 1. Introduction

Superconductors and superfluids display a wide range of critical temperatures, gap structures and transport behaviors, yet many empirical regularities survive across families. Classical Type-I and Type-II superconductors, high-pressure hydrides, oxide and Fe-based compounds, and superfluid helium all exhibit sharp transitions and non-trivial scaling laws that are not easily reduced to a single microscopic mechanism. A natural question is whether there exists a *coarse-grained* description that can organize these systems into a small number of classes, using only a few parameters and without per-material tuning of exponents.

The DOFT (Differential Oscillator Field Theory) framework proposes such a coarse-grained description. In DOFT, clusters and “grainy” structures are encoded through resonance locks between oscillators, leading to a discrete set of admissible exponents. Rather than fitting continuous exponents freely, the theory prescribes a small library of integer and rational exponents arising from normal locking conditions. In this sense, DOFT is not used here as a microscopic theory of superconductivity, but as a **locking grammar** that constrains how effective fingerprints are built from a limited set of exponents.

In this work we take a deliberately conservative stance: we do not attempt to prove DOFT as a fundamental theory. Instead, we ask a narrower and more practical question: **given** a DOFT-inspired locking grammar and a minimal correction law, can we describe a heterogeneous collection of superconductors and superfluid helium with a small, fixed set of parameters? Concretely, each material is mapped to either an integer fingerprint or a rational fingerprint \(q\), whose exponents are fixed by the locking rules and *not* adjusted per material. A universal correction law with parameters \(\Gamma\) and \(\eta\) is then calibrated on classical metals under weak constraints (\(\Gamma \ge 0\), \(\eta \ge 0\)) and reused across families.

Using this setup, we perform a series of statistical tests on published data for classical superconductors, high-pressure hydrides, oxides, Fe-based superconductors and superfluid helium. We first show that the calibration on classical metals yields \(\Gamma \approx 0\) and a robust positive \(\eta\), as confirmed by bootstrap and leave-one-out influence analysis. We then examine the stability of integer fingerprints across runs and their ability (or lack thereof) to distinguish Type-I from Type-II superconductors. Non-parametric tests (Kruskal–Wallis, Mann–Whitney, Kolmogorov–Smirnov) and multiple-comparison corrections are used to quantify differences between families in specific exponents and in the rational fingerprint \(q\). Finally, we explore a sub-network protocol for detecting cluster-level couplings, and we assess residual diagnostics by family to test whether the proposed correction law removes systematic drift.

Our results indicate that (i) a single dissipative scale \(\eta\) suffices to describe classical metals, with \(\Gamma\) effectively null; (ii) Type-I and Type-II superconductors share a common integer fingerprint within the noise of our data; (iii) high-pressure hydrides separate from classical metals along a specific integer exponent, yielding a moderate but statistically robust effect; (iv) superfluid helium and high-pressure hydrides occupy distinct, tightly clustered regions in rational \(q\)-space, clearly separated from oxides and Fe-based superconductors; and (v) within the current protocol, cluster-level couplings are not required to explain the data, except for a small but detectable effect in MgB₂. We interpret these findings as evidence that a DOFT-inspired locking grammar, combined with a minimal correction law, can serve as a compact organizational scheme for diverse superconducting and superfluid systems, while remaining fully testable and falsifiable as new data become available.

## 2. Methods

### 2.1. Data and families

We compile published data for superconductors and superfluid helium into a canonical dataset (`materials_master_v5.csv`). Each row corresponds to one material and includes, at minimum, its name, family label (classical metal, high-pressure hydride, oxide, Fe-based, superfluid), critical temperature \(T_c\), and additional descriptors as needed. Families are defined phenomenologically rather than by microscopic mechanism, in order to test whether a single DOFT-inspired grammar can span heterogeneous systems.

Classical metals (both Type-I and Type-II) are used as the calibration set for the correction law parameters \(\Gamma\) and \(\eta\). High-pressure hydrides, oxides, Fe-based superconductors and superfluid helium are then treated as *test families* for the fingerprints, with no additional per-material tuning of exponents.

### 2.2. DOFT-inspired fingerprints

In the DOFT framework, resonance locks between oscillators restrict the allowed exponents to a discrete set of integer and rational values. We use this locking grammar to define two types of fingerprints:

- **Integer fingerprints**, characterized by a small set of fixed exponents \((\alpha_1, \alpha_2, \dots)\) that encode “integer-locking” behavior. These are used primarily for classical metals and high-pressure hydrides.
- **Rational fingerprints**, characterized by a rational parameter \(q\) and associated exponents of the form \(p/q\) with small denominator \(q\). These are used for softer or hybrid systems such as oxides, Fe-based superconductors and superfluid helium.

Crucially, the exponents themselves are **pre-specified** by the DOFT locking rules and not fitted separately for each material. The fingerprints are therefore highly constrained: given a family and the locking grammar, the exponents are fixed and only a small number of continuous parameters remain to be calibrated globally.

### 2.3. Universal correction law

We introduce a universal correction law with parameters \(\Gamma\) and \(\eta\) to account for curvature and dissipative memory effects, respectively. For each material, the observable of interest (e.g. a normalized transition or scaling quantity) is modeled as a function of the fingerprint exponents and the layer index or effective distance \(d\), modified by a correction term

\[
\Delta(d) \approx \Gamma f_\Gamma(d) + \eta f_\eta(d),
\]

where \(f_\Gamma\) and \(f_\eta\) are fixed functions determined by the DOFT-inspired ansatz. The calibration is performed on classical metals only, under weak inequality constraints

\[
\Gamma \ge 0, \quad \eta \ge 0.
\]

We use constrained least-squares to obtain point estimates \(\hat{\Gamma}, \hat{\eta}\), and we quantify uncertainty and robustness through bootstrap and leave-one-out analyses (Section 2.4). Once calibrated, \(\hat{\Gamma}\) and \(\hat{\eta}\) are held fixed for all other families.

### 2.4. Bootstrap and influence analysis

To assess the stability of the calibration on classical metals, we perform a non-parametric bootstrap with \(N = 500\) resamples. For each bootstrap sample we re-fit the correction law, obtaining distributions for \(\Gamma\) and \(\eta\). We report percentile-based confidence intervals and examine the skewness and width of these distributions.

We also carry out a leave-one-out (LOO) influence analysis: each classical metal is removed in turn, the model is re-fitted, and the resulting \(\Gamma\) and \(\eta\) estimates are compared to the full-sample values. This allows us to detect whether any single material exerts disproportionate influence on the calibration.

### 2.5. Statistical tests for family separation

We use standard non-parametric tests to quantify differences between families in their fingerprint parameters:

- **Kruskal–Wallis (KW)** tests for global differences in specific exponents (e.g. `exp_d_7`, `exp_a_2`) across multiple families.
- **Mann–Whitney U** tests for pairwise comparisons, in particular between high-pressure hydrides and classical metals along selected integer exponents.
- **Kolmogorov–Smirnov (KS)** tests for differences in the distribution of rational \(q\) between families (e.g. superfluid helium vs oxides and Fe-based superconductors).

To control for multiple comparisons when scanning over several exponents, we apply false discovery rate (FDR) correction (Benjamini–Hochberg). Effect sizes are reported alongside p-values, using Cliff’s delta for two-sample tests where appropriate.

### 2.6. Sub-network protocol and cluster coupling

To probe possible cluster-level couplings, we define a sub-network protocol based on a contrast \(C_{AB}\) between two channels (e.g. \(\sigma\) and \(\pi\) bands in MgB₂). A coupling parameter \(\kappa \ge 0\) is introduced only when the contrast suggests potential mixing. The protocol is intentionally conservative: \(\kappa\) is included as an extra parameter only if it reduces residual drift in a statistically meaningful way; otherwise, \(\kappa\) is set to zero.

This procedure is applied across the dataset, with particular attention to systems such as MgB₂ and FeSe where multi-band or multi-channel behavior is expected. We treat the detection of \(\kappa \approx 0\) as informative **negative evidence** for cluster-level coupling under the present DOFT-inspired ansatz.

### 2.7. Residual diagnostics

Finally, we analyze residuals by family after applying the universal correction law and the fingerprint exponents. For each family we examine:

- the mean and variance of residuals;
- residuals as a function of layer distance or effective coordinate \(d\);
- the presence or absence of systematic drift or structure in residual plots.

The goal is to test whether the proposed parametrization removes coherent trends, leaving residuals centered around zero with reasonable dispersion. Clean residuals, especially in high-pressure hydrides, are interpreted as support for the sufficiency of the DOFT-inspired grammar and the universal correction law.

## 3. Results

### 3.1. Calibration on classical metals: Γ ≈ 0, η > 0

Applying the constrained calibration procedure to classical metals yields

- \(\hat{\Gamma} \approx O(10^{-17})\), effectively indistinguishable from zero within numerical precision and bootstrap uncertainty;
- \(\hat{\eta} \approx 4 \times 10^{-5}\), a robust positive value.

Bootstrap resamples (N = 500) show a narrow, unimodal distribution for \(\eta\) with confidence intervals that remain strictly positive. The distribution for \(\Gamma\) is sharply peaked around zero and consistent with the non-negativity constraint acting as a hard bound rather than truncating a broader bulk. Leave-one-out influence analysis confirms that no single classical metal dominates the calibration: removing any individual material shifts \(\eta\) by at most a small fraction of its estimated value, and \(\Gamma\) remains negligible in all cases.

These results support the interpretation that curvature in the correction law is not required for classical metals (\(\Gamma \approx 0\)), and that a single dissipative scale \(\eta\) suffices to describe this family under the DOFT-inspired ansatz.

### 3.2. Stable integer fingerprints and Type-I vs Type-II comparison

Using the integer-locking grammar, we assign fixed exponent sets to classical superconductors, yielding stable fingerprints for Type-I and Type-II families. Across repeated runs and resampling, the mean exponents remain close to

- Type-I: approximately \((1.55, 0.82, 0.52, 0.41)\),
- Type-II: approximately \((1.93, 0.59, 0.52, 0.41)\),

with small run-to-run variation.

A global Kruskal–Wallis test across multiple families detects significant differences in at least one exponent (notably `exp_d_7` and, in the latest runs, `exp_a_2`), indicating that the integer fingerprints carry meaningful structure at the family level. However, post-hoc pairwise comparisons between Type-I and Type-II superconductors do *not* yield statistically significant separation once multiple comparisons are taken into account. Within the uncertainty of our data, Type-I and Type-II share a common integer fingerprint class.

This supports the idea that, in the DOFT-inspired picture, classical superconductors belong to a single integer-locking family, despite their conventional BCS classification differences.

### 3.3. High-pressure hydrides vs classical metals: a fingerprint signature

In contrast to the Type-I/Type-II comparison, high-pressure hydrides (HP) exhibit a clear separation from classical metals along a specific integer exponent, denoted `exp_d_7`. A Mann–Whitney U test on this exponent yields an FDR-corrected p-value of order \(10^{-4}\), with a moderate effect size (Cliff’s delta ≈ −0.30). The sign and magnitude of the effect indicate a systematic shift in the HP fingerprint relative to classical metals, rather than a difference driven by a few outliers.

We interpret this as a **quantitative fingerprint signature** for the high-pressure class within the integer-locking space. The same DOFT-inspired grammar that fails to distinguish Type-I from Type-II superconductors is thus able to separate HP from classical metals along a constrained axis, with no additional exponents or per-material fitting.

### 3.4. Rational fingerprints: superfluid helium and high-pressure hydrides

For softer or hybrid systems we consider rational fingerprints encoded by a parameter \(q\). Superfluid helium occupies a distinctly low-\(q\) regime, with values typically in the range \(q \approx 2\)–2.3. Oxides and Fe-based superconductors populate higher-\(q\) regions, and Kolmogorov–Smirnov tests reveal highly significant differences between the superfluid and these families (p ≪ 0.01).

High-pressure hydrides, when described in the rational fingerprint space, cluster around \(q \approx 5.85\) with tight confidence intervals. The narrow spread suggests an internally coherent DOFT class: HP systems share a common rational-locking pattern that is distinct both from superfluid helium and from oxide/Fe-based families.

Overall, the rational fingerprints support a picture in which superfluid helium and high-pressure hydrides occupy separate, tightly clustered regions in \(q\)-space, clearly separated from oxides and Fe-based superconductors. This is consistent with the broader DOFT claim that soft/hybrid families are characterized by rational \(p/q\) locking with small denominators.

### 3.5. Sub-network protocol and cluster coupling

Applying the sub-network protocol across the dataset, we find that the cluster-coupling parameter \(\kappa\) is statistically compatible with zero for all systems except MgB₂. In MgB₂, the contrast \(C_{AB}\) between the two channels (e.g. \(\sigma\) and \(\pi\)) is substantial (\(C_{AB} \approx 1.59\)), and allowing a small but non-zero \(\kappa \approx 6 \times 10^{-3}\) reduces residual drift. For FeSe and other candidate multi-band systems, the inferred \(\kappa\) remains negligible within uncertainty.

This pattern provides negative but useful evidence: under the current DOFT-inspired rules and correction law, no significant cluster-level couplings are required to explain the data, except for MgB₂ as a mixed-channel outlier. The detection of a small, localized coupling in MgB₂ is consistent with its known multi-band character.

### 3.6. Residual diagnostics by family

Residual diagnostics by family show means close to zero and no explosive dispersion after applying the universal correction law and the family-specific fingerprints. High-pressure hydrides in particular exhibit clean residuals, with little visible structure when plotted as a function of layer distance or effective coordinate \(d\).

The absence of systematic drift in residuals suggests that the combination of:

1. a universal dissipative scale \(\eta\),
2. a null or negligible curvature parameter \(\Gamma\), and
3. the discrete locking grammar for exponents,

is sufficient to capture the dominant structure in the data for the systems considered here. Additional ad-hoc parameters or exponents are not needed within the current level of noise and sampling.

## 4. Discussion

### 4.1. DOFT-inspired grammar as a phenomenological model

The central question of this study is not whether DOFT is a complete microscopic theory of superconductivity or superfluidity, but whether its locking grammar can serve as a useful phenomenological parametrization. By constraining exponents to a discrete set derived from resonance locks, and by calibrating only a small number of continuous parameters globally, we obtain a compact description that remains falsifiable: if the data were incompatible with the grammar, this would show up in poor fits, unstable calibrations or structured residuals.

Our findings indicate that the DOFT-inspired grammar clears several non-trivial tests. Classical metals can be described with a single dissipative scale \(\eta\) and negligible curvature \(\Gamma\); Type-I and Type-II superconductors share a common integer fingerprint; high-pressure hydrides separate from classical metals along a specific constrained exponent; superfluid helium and high-pressure systems occupy distinct, coherent regions in rational \(q\)-space; and cluster-level couplings are not required except for a localized effect in MgB₂. None of these patterns are guaranteed a priori, given the tight constraints on exponents and parameters.

### 4.2. Comparison to simpler baselines

A natural concern is whether the observed patterns could be captured equally well, or better, by simpler models without a discrete locking grammar—for example, by fitting arbitrary power laws or flexible scaling forms to each family. From a statistical perspective, such baselines typically introduce more degrees of freedom, allowing per-family or even per-material exponents, at the cost of reduced falsifiability.

In contrast, the DOFT-inspired approach explicitly restricts the exponent space and emphasizes cross-family reuse of parameters. The fact that family-level differences (e.g. between high-pressure hydrides and classical metals) appear along specific constrained exponents, while other expected distinctions (Type-I vs Type-II) do not, is consistent with a picture in which the grammar is neither overly rigid nor overly flexible. A more detailed model comparison against specific baselines—quantifying goodness-of-fit, parameter count and predictive performance—remains an important direction for future work.

### 4.3. Limitations

Several limitations should be kept in mind. First, the number of materials per family, especially for high-pressure hydrides and some oxide or Fe-based subsets, is still modest. This limits the statistical power of some tests and increases sensitivity to selection effects in the available literature. Second, the present implementation relies on a particular choice of observables and derived quantities; alternative constructions might interact differently with the DOFT-inspired grammar. Third, the sub-network protocol for detecting cluster-level couplings is intentionally conservative and may miss more subtle forms of mixing.

These limitations do not invalidate the patterns reported here but constrain how strongly they should be interpreted. We view the present work as an initial test of a structured parametrization, rather than as a definitive classification of all superconductors and superfluids.

### 4.4. Outlook

Looking forward, several extensions suggest themselves. On the data side, expanding the dataset with additional high-pressure hydrides, unconventional superconductors and more detailed measurements of superfluid helium would sharpen the statistical tests and better probe the boundaries of the locking grammar. On the modeling side, systematic comparisons with simpler baselines and with alternative locking schemes would clarify which aspects of the results are specific to DOFT and which are more generic.

More broadly, the approach taken here—using a discrete locking grammar and a small number of global parameters to organize diverse systems—may be applicable beyond superconductivity and superfluidity. Any domain where layered or clustered structures emerge from oscillator-like dynamics could provide a testing ground for DOFT-inspired fingerprints. Whether or not DOFT ultimately survives such tests, the methodology of tightly constrained, falsifiable parametrizations may be a useful addition to the toolbox for organizing complex condensed-matter phenomena.

## 5. Conclusions

We have tested a DOFT-inspired fingerprint model on a heterogeneous collection of superconductors and superfluid helium, using a tightly constrained set of integer and rational exponents combined with a universal correction law. The main conclusions are:

1. A single dissipative parameter \(\eta \approx 4 \times 10^{-5}\) suffices to describe classical metals, while the curvature parameter \(\Gamma\) is effectively null.
2. Integer fingerprints are stable across runs and do not significantly distinguish Type-I from Type-II superconductors, indicating a shared integer-locking class.
3. High-pressure hydrides separate from classical metals along a specific constrained exponent (`exp_d_7`), yielding a statistically robust, moderate effect that serves as a quantitative fingerprint signature.
4. In rational \(q\)-space, superfluid helium and high-pressure hydrides occupy distinct, tightly clustered regions that are clearly separated from oxides and Fe-based superconductors.
5. Cluster-level couplings, as encoded by \(\kappa\), are not required to explain the data under the present ansatz, except for a small but detectable effect in MgB₂.
6. Residual diagnostics by family show means close to zero and reasonable dispersion, especially for high-pressure hydrides, suggesting that no additional ad-hoc parameters are needed at the current level of detail.

Taken together, these results support the view that a small DOFT-inspired parameter set, combined with a discrete locking grammar, can provide a compact and testable organizational scheme for superconducting and superfluid systems. The framework remains falsifiable as new data and more stringent baselines are brought to bear, and it invites further exploration of resonance-based fingerprints in other complex materials.

## 6. Data and code availability

All data and code used in this study are intended to be made publicly available in this repository. The canonical input dataset will be provided as `data/raw/materials_master_v5.csv`, and the main analysis pipeline will be implemented in the `src/` directory. Processed results (calibration outputs, fingerprints, statistical test summaries) and the figures and tables corresponding to the manuscript will be stored under `data/processed/` and `results/` respectively.

Once finalized, this draft can be updated with a permanent DOI (e.g. via Zenodo or OSF) and a formal citation entry.

## Acknowledgements

The author thanks informal collaborators and the broader condensed-matter community for maintaining open data and discussions that made this exploratory study possible. Any errors or misinterpretations remain the responsibility of the author.
