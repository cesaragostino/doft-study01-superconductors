# References – DOFT Study 01 (Superconductors & Superfluid Helium)

This bibliography collects the core sources we draw on for DOFT Study 01:
delay/memory formalisms, superconductivity and superfluid helium, high-pressure
hydrides and unconventional superconductors, signal models, and reproducible
computational practice.

We use bracketed tags like **[DDE-1]**, **[SC-2]**, **[REP-1]** in the text.  
Each tag maps to an entry below.

---

## A. Delays, memory, and coarse-graining

These works motivate the use of delay/memory kernels and coarse-grained
descriptions, which underlie the DOFT-inspired correction law.

**[DDE-1]** Jack K. Hale & S. M. Verduyn Lunel, *Introduction to Functional Differential Equations*, Springer, 1993.  
**[DDE-2]** Richard Bellman & Kenneth L. Cooke, *Differential-Difference Equations*, Academic Press, 1963.  
**[MEM-1]** Ryogo Kubo, “The fluctuation–dissipation theorem,” *Rep. Prog. Phys.* **29** (1966) 255–284.  
**[MEM-2]** Hazime Mori, “Transport, collective motion, and Brownian motion,” *Prog. Theor. Phys.* **33** (1965) 423–455.  
**[MEM-3]** Robert Zwanzig, *Nonequilibrium Statistical Mechanics*, Oxford Univ. Press, 2001.  

---

## B. Superconductivity and superfluid helium

These references provide standard background on superconductivity, superfluid
helium, and related effective descriptions.

**[SC-1]** Michael Tinkham, *Introduction to Superconductivity*, 2nd ed., McGraw-Hill, 1996.  
**[SC-2]** J. Bardeen, L. N. Cooper & J. R. Schrieffer, “Theory of Superconductivity,” *Phys. Rev.* **108** (1957) 1175–1204.  
**[SC-3]** A. J. Leggett, *Quantum Liquids: Bose Condensation and Cooper Pairing in Condensed-Matter Systems*, Oxford Univ. Press, 2006.  
**[SC-4]** D. R. Tilley & J. Tilley, *Superfluidity and Superconductivity*, 3rd ed., Hilger, 1990.  

---

## C. High-pressure hydrides, oxides, and Fe-based superconductors

These works give context for high-\(T_c\) hydrides and other unconventional
superconductors considered in the dataset.

**[HP-1]** N. W. Ashcroft, “Hydrogen Dominant Metallic Alloys: High Temperature Superconductors?,” *Phys. Rev. Lett.* **92** (2004) 187002.  
**[HP-2]** Russell J. Hemley & H. K. Mao, “Progress in high-temperature superconductivity in hydrides at high pressure,” *Rev. Mod. Phys.* (review articles / overviews as appropriate).  
**[OX-1]** J. G. Bednorz & K. A. Müller, “Possible high-\(T_c\) superconductivity in the Ba–La–Cu–O system,” *Z. Phys. B* **64** (1986) 189–193.  
**[FE-1]** Y. Kamihara *et al.*, “Iron-based layered superconductor La[O\(_{1-x}\)F\(_x\)]FeAs (x = 0.05–0.12) with \(T_c = 26\) K,” *J. Am. Chem. Soc.* **130** (2008) 3296–3297.  

*(You can refine HP-2 to the specific review you actually use.)*

---

## D. DOFT and locking grammar (context)

These references are specific to the DOFT program and the locking grammar used
to constrain exponents in this study.

**[DOFT-1]** C. Agostino, “Delayed Oscillator Field Theory (DOFT): Mother Frequency and Thermal Memory Shift,” internal notes / preprint (in preparation).  
**[DOFT-2]** C. Agostino, “DOFT Manifesto v1.8,” project document describing the broader delayed-oscillator program and locking rules (2025).  

*(Update [DOFT-1]/[DOFT-2] once there is a public preprint or DOI.)*

---

## E. Chaos, self-organization, and criticality

These works provide background on self-organized criticality and lack of
self-averaging in disordered systems, relevant for interpreting scaling and
cluster formation.

**[SOC-1]** Per Bak, Chao Tang & Kurt Wiesenfeld, “Self-organized criticality,” *Phys. Rev. A* **38** (1988) 364–374.  
**[NAV-1]** A. Aharony & A. B. Harris, “Absence of self-averaging and universal fluctuations in random systems near critical points,” *Phys. Rev. Lett.* **77** (1996) 3700–3703.  

---

## F. Signal models, exponential sums, and uncertainty

These underpin aspects of the fingerprint construction and the treatment of
exponential components.

**[PRN-1]** G. R. de Prony, “Essai expérimental et analytique…,” *J. École Polytech.* **1** (1795) 24–76.  
**[DSP-1]** Alan V. Oppenheim & Ronald W. Schafer, *Discrete-Time Signal Processing*, 3rd ed., Pearson, 2010.  
**[PRN-2]** M. R. Osborne & G. K. Smyth, “A modified Prony algorithm for exponential function fitting,” *SIAM J. Sci. Comput.* **16** (1995) 119–138.  

---

## G. XY model, Anderson–Higgs, and gauge/phase context

These are optional but useful for interpreting phase-locking, gauge-like
structures, and coarse-grained order parameters.

**[XY-1]** B. D. Josephson, “Possible new effects in superconductive tunnelling,” *Phys. Lett.* **1** (1962) 251–253.  
**[AH-1]** P. W. Anderson, “Plasmons, gauge invariance, and mass,” *Phys. Rev.* **130** (1963) 439–442.  
**[RG-1]** K. G. Wilson, “Renormalization group and critical phenomena,” *Phys. Rev. B* **4** (1971) 3174–3183.  

---

## H. Methodological and reproducibility references

These references motivate the way the repository is structured and how we handle
computational reproducibility.

**[REP-1]** Victoria Stodden, “Reproducible Research: Tools and Strategies for Scientific Computing,” *Comput. Sci. Eng.* **11** (2009) 11–12.  
**[REP-2]** Sandve et al., “Ten simple rules for reproducible computational research,” *PLoS Comput. Biol.* **10** (2014) e1003285.  

---

*Note:* This file is specific to DOFT Study 01. For a broader list of references
covering the full DOFT program (including emergent gravity, Rydberg/QDT, etc.),
see the main DOFT repository.
