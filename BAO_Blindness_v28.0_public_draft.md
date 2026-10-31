# Per-Probe Decomposition of Bimetric Cosmology: BAO Blindness and a Baseline Scalar-Sector Diagnostic
Aleksandr Chernikov
*Independent Researcher*
1 November 2026 (v28.0, public working draft)
## Abstract
We present a conditional, background-level, per-probe analysis of the Hassan–Rosen bimetric benchmark used by Högås & Mörtsell (2025), using DESI DR2 BAO data (13 points), Pantheon+ SNe Ia (1588 light curves after the z > 0.01 cut), and the compressed Planck 2018 + ACT DR6 likelihood. At the baseline point, the BAO difference from the fitted ΛCDM reference is small (Δχ² = 0.188 with the full covariance). We use “BAO blindness” only in this restricted sense of statistical indistinguishability for the tested baseline and grids, not as absence of physical coupling or a property of the entire bimetric theory. The mechanism is the small late-time bimetric correction to E²(z), together with profiling of the scale f = c/(H₀r_d); profiling is minimisation over the nuisance parameter, not Bayesian marginalisation.

For Pantheon+, the corresponding baseline difference is Δχ² = 6.50 with the full STAT+SYS covariance and 15.57 with the diagonal covariance, the latter being an upper-bound diagnostic. These are per-probe comparisons at the fixed benchmark, not a joint model-selection significance, evidence, or Bayes factor. Across the reported 20-point (α, β) grid and the stated local analytic approximation, the predicted Ω_m is nearly rigid; this restricted result is not a proof over the full parameter space. The compressed-CMB geometric diagnostic gives a small baseline penalty and a large penalty as θ increases to 30° in the P1 convention (Δχ² = 0.04097 → 2369.29001 against the free three-parameter reference). Because that likelihood is GR-calibrated and no bimetric Boltzmann/ISW calculation is included, the physical error of transferring it to bimetric gravity remains unquantified.

At the same baseline and on the adopted finite branch, the linear/WKB scalar-sector diagnostic locates a gradient-instability boundary z* = 7.880536 and gives a large accumulated growth exponent. This is a baseline perturbative result, not a viability map for the full model class. A subsequent joint P3 analysis with exact Fermi–Dirac neutrinos is under methodological audit; no accepted P3 posterior or model-comparison result is used in this paper.
**Keywords: bimetric gravity, DESI BAO, dark energy, modified gravity, CMB, supernovae**

## Public-draft status

This is a public working draft released on 1 November 2026 while the methodological audit continues. It supersedes manuscript v21.6 as the current public account of the project, but it is not a final model-selection paper.

- The reported numerical results belong to the stated legacy per-probe conventions and scanned domains. Stronger joint, Bayesian, or full-CMB interpretations are outside their evidence.
- The developing P3 target uses exact Fermi–Dirac neutrinos and a restricted computational domain. It has no accepted production MCMC result here and does not support the numerical claims of this draft.
- Profiling of f and M is not full marginalisation over those nuisance parameters. The quantitative effect of that statistical choice on scientific conclusions remains under audit.
- Reproducibility is version-specific. The executable repository archive v10.6 accompanies manuscripts v21.5/v21.6 and does not reproduce every updated number in v28.0. The documentary code archive v13.3 is included in this public package; the distinction is specified in `REPRODUCIBILITY.md`.
## 1. Introduction
Bimetric gravity (Hassan & Rosen 2012 [Hassan2012]) extends general relativity by a dynamical second metric f_μν coupled to the physical metric g_μν through a Fierz–Pauli mass term, while remaining ghost-free at the fully nonlinear level; it admits viable cosmological solutions on an FLRW background (von Strauss et al. 2012 [vonStrauss2012]; Comelli et al. 2012 [Comelli2012]).

This work extends the Högås & Mörtsell (2025) [Hogas2025] setup by decomposing the constraining power of each cosmological probe individually. The decomposition delivers a diagnosis that the joint analysis alone cannot: BAO is statistically indistinguishable from ΛCDM at the baseline parameters (Δχ² = 0.188, full covariance), while SNe Ia and the compressed CMB geometric proxy are sensitive. We therefore propose BAO blindness — in the precise sense of statistical indistinguishability (Δχ² ≪ 3.84), not absence of coupling — as a diagnostic tool for the bimetric class: which probes constrain the parameters, which are blind, and why.

Why bimetric gravity: it is a consistent, ghost-free extension of GR with a rigid prediction structure — the quartic constraint predicts Ω_m from the model parameters rather than fitting it. The sensitivity/blindness pattern of the three probes is therefore a testable diagnostic hypothesis for the tested model settings, not an established property of the full class.

Högås & Mörtsell (2021) [Hogas2021] derived the B-parameter formalism (observational constraints in the companion paper [Hogas2021b]), and Högås & Mörtsell (2025) provided the full Friedmann equation, Ω_DE(z), and quartic solver for the scale factor ratio y(z). Their analysis combined BAO, CMB, and SNe Ia data through MCMC to obtain joint constraints on the bimetric parameters, but did not decompose the constraining power of individual probes. Table 1 summarises the key differences. Our original contributions are: (1) the BAO χ² pipeline with profiling of the scale f = c/(H₀r_d); (2) numerical B-parameter values; (3) the per-probe Δχ² decomposition revealing BAO blindness and SNe/CMB sensitivity; (4) the decomposition of the BAO Δχ² into the Ω_m-mismatch and pure E²(z)-correction contributions; (5) the round-trip consistency test (10⁻¹³ on B-parameters, 10⁻¹⁵ on E²(0)); (6) the CMB-only θ sensitivity as an isolated result (0.04 → 2369, P1 convention); (7) the full-covariance SNe assessment (Δχ² = 6.50), which closes the leading quantitative limitation of the diagonal analysis; and (8) a quantitative scalar-sector stability analysis at the baseline point: the gradient-instability boundary z* and the accumulated growth factor Γ(k) (§6).

Lüben, Schmidt-May & Smirnov (2020) [Luben2020] demonstrated that the early universe in bimetric theory is screened by the Vainshtein mechanism on an FLRW background, with the spin-2 mass setting the cosmological Vainshtein scale. In a global analysis using CMB, Cepheid, SNe Ia, and BAO data, they showed that a minimal β_0β_1β_4 bimetric model can resolve the H₀ tension via an effective phantom dark energy component (w_DE < -1) in the redshift interval 1 ≲ z ≲ 10. Their finding that current BAO data do not distinguish bimetric models from ΛCDM is consistent with and complementary to our analysis: we quantify this indistinguishability per probe and identify its origin in the smallness of the bimetric correction at late times.

The restricted rigidity result underlying this part of the paper: the predicted Ω_m = 0.304 is independent of the potential parameters to |∂Ω_m/∂α| ≈ 4 × 10⁻⁹ and |∂Ω_m/∂β| ≈ 10⁻¹² across the scanned grid — the Ω_m tension is structural only in the restricted grid-and-track sense defined in §7.5.

```
Table 1. Comparison with Högås & Mörtsell (2025).

Aspect                 H&M 2025               This work
---------------------  ---------------------  ------------------------------------
Data                   BAO + CMB + SNe Ia     Each probe individually
                       (combined)
Method                 MCMC (joint)           Grid scan + Brent (per probe)
Ω_m                    Derived from quartic   Predicted by quartic
                       (varies in MCMC)
B-parameters           Not reported           Numerical (§2.4)
Δχ² decomposition      Not performed          Per probe (Table 2)
BAO blindness          Not quantified         Δχ² = 0.188 (full, incl. Ω_m
                                              penalty; 0.021 pure E²(z)
                                              correction, see §3.2), flat across
                                              20 (α, β) combinations and 4
                                              decades in m_FP/H₀
SNe sensitivity        Not quantified         Δχ² = 15.57 (diag), 6.50 (full cov);
                                              1.49/0.96 (pure); 11.58 with
                                              σ_int = 0.12
CMB sensitivity        Not isolated           absolute χ² = 0.041–2369.29; P1 Δχ² = 0.041–2369.29
to θ
θ sensitivity          Not provided           flat (0.204→0.198) in the
(BAO only)                                    self-consistent scan; 56.1 at
                                              θ = 31.5° at fixed Ω_m
                                              (block-diagonal covariance)
Round-trip test        Not performed          10⁻¹³ (B-param), 10⁻¹⁵ (E²)
```

## 2. Theoretical Framework
### 2.1. Bimetric Friedmann Equation
The bimetric Friedmann equation for the physical metric g_μν takes the form (Högås & Mörtsell 2025, Eq. A12): E²(z) ≡ H²(z)/H₀² = Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_DE(z), where Ω_DE(z) = Ω_Λ - sin²θ · m_FP² · (1-y)[1 + α(1-y) + β/3(1-y)²], where θ is the mixing angle, m_FP is the Fierz–Pauli mass, α and β are potential parameters, and Ω_Λ is the effective cosmological constant, treated as a free model parameter (following Högås & Mörtsell 2025, who parameterize it via ω_Λ = Ω_Λ h²); throughout this work Ω_Λ plays the role of the joint-fit anchor of the three probes (see §7.5, where we show that moving this anchor along the closure degeneracy ray costs Δχ² ≈ 47.9 in total). At the baseline parameters (α = 0, β = 3, m_FP/H₀ = 25, θ = 6.1°), Ω_Λ ≈ 0.699. In the GR limit (y → 1, θ → 0), Ω_DE → Ω_Λ and the standard Friedmann equation is recovered. The physical interpretation is that Ω_DE receives a negative contribution proportional to sin²θ · m_FP², representing the energy density of the massive graviton. As y deviates from 1 (the GR limit), this correction grows, modifying the expansion history.

*Note on the late-time effective form.* Expanding Ω_DE(z) about the de Sitter attractor y → 1 yields the effective law E²(z) ≈ Ω_Λ + cos²θ·Ω_m(z), in which the entire late-time bimetric modification enters through the single rescaling Ω_m(z) → cos²θ·Ω_m(z). This law is not new to this work: it is the exponential-expansion form of Högås & Mörtsell (2021) [Hogas2021] (their dark-degeneracy relation, with the identifications Ω_m0 = (1 − ω_leff)/cos²θ and Ω_de = ω_leff − sin²θ·Ω_m). Our independent derivation — a series expansion in ε = 1 − y with resummation — rediscovered the same result; we therefore regard the agreement as a validation of our expansion method, not as an original contribution. What remains specific to this work is the quantitative redshift map of the approximation's accuracy: at the baseline point the residual |E²_bim − E²_eff|/E²_bim reaches only 0.24% at z = 6.4 and 0.64–0.75% at z = 11–12 when the radiation term is omitted from the effective law (0.019% at z = 6.4 with the Ω_r term retained, the residual in this convention passing through zero at z = 4.74), crosses 1% only at z ≈ 14.5 (z ≈ 31 with Ω_r retained), and asymptotes to the matter-era plateau sin²θ ≈ 1.13% from below — there is no breakdown at z ≈ 11–12, and the effective law remains a sub-percent approximation well past the Vainshtein boundary at z ≈ 12 (§7.1, Table A2). This map is confirmed by four independent blind external re-implementations (verbatim runs, Appendix C).
### 2.2. Quartic Equation for y(z)
The scale factor ratio y(z)[^1] satisfies a quartic equation (Högås & Mörtsell 2025, A13; corresponding to Eq. 5.5 of Högås & Mörtsell 2021 [Hogas2021]): a_0 + a_1 y + a_2 y² + a_3 y³ + a_4 y⁴ = 0, with coefficients per Eq. A13 of Högås & Mörtsell (2025); written explicitly through the B-parameters of §2.4 and t ≡ tan²θ (H0 = 1 units): a_4 = t·B_3, a_3 = 3t·B_2 − B_4, a_2 = 3t·B_1 − 3B_3, a_1 = t·(B_0 + 3Ω_tot(z)) − 3B_2, a_0 = −B_1, where Ω_tot(z) = Ω_m(1+z)³ + Ω_r(1+z)⁴ and m_FP² is in units of H₀². The coefficients follow the B-parameter formalism of Högås & Mörtsell (2021, 2025), solving for y = a_f/a_g directly. We solve numerically using Brent's method on the interval [10⁻⁶, 1.0] (the lower end chosen to avoid the singular limit while containing the full physical branch), which robustly finds the physical (finite) branch: y → 1 as z → 0 (GR limit) and y → 0 as z → ∞ (early universe), consistent with the finite branch definition of Högås & Mörtsell (2021, Sec. 5.1) [Hogas2021]. The infinite branch, where y grows with z, is excluded due to the Higuchi ghost (Högås & Mörtsell 2021, Sec. 5.1). The solution is verified by computing E²(0) = 1, with a residual of ~10⁻¹⁵ (2 × 10⁻¹⁵ in the canonical run; environment-dependent), consistent with the self-consistency condition of Högås & Mörtsell 2021 [Hogas2021], Eq. 5.6. Note that E²(0) = 1 exactly by the self-consistency condition Ω_m = 1 - Ω_r - Ω_DE(0); the residual ~10⁻¹⁵ reflects the numerical precision of Brent's method in solving the quartic equation and serves as a numerical consistency check.
*Figure 1: Fractional bimetric correction to the expansion history, δE² = E²_bim/E²_ΛCDM − 1*
### 2.3. Self-Consistent Ω_m
The bimetric model predicts Ω_m from the quartic constraint, rather than fitting it as a free parameter. At the baseline parameters (α = 0, β = 3, m_FP/H₀ = 25, θ = 6.1°), the quartic yields y(0) = 0.99951 and Ω_m = 0.30435 (BAO/SNe-pipeline convention for Ω_r and closure; see the note on the Ω_m conventions in Appendix A). The fractional bimetric correction to E²(z) is shown in Fig. 1 and tabulated in Table A2.
**Analytic closure remark (α = 0, β = 3).** The approximate constancy of C·Ω_m, where C = cos²θ, is a finite-parameter analogue of the bimetric realization of Kunz's dark degeneracy [Kunz2009]; the exact screened-limit relation was derived by Högås & Mörtsell (2021) [Hogas2021]. We quantify the controlled departure from that limit rather than claim a new degeneracy or symmetry. With S = sin²θ, L = Ω_Λ, e = 1 − y(0), Ω_t0 = Ω_m + Ω_r and M = m²(e + e³), the exact z = 0 closure is

  Ω_m = 1 − L − Ω_r − Ω_ν + S·M,
  C·Ω_m = Ω_0 + S·(C·M − Ω_0),   Ω_0 = 1 − L − Ω_r − Ω_ν.

Here Ω_m is the nonrelativistic-matter density excluding neutrinos, as in the CMB convention. Defining D₁ = m² + Ω_t0 − 2L, the quartic gives e₁ = Ω_t0/D₁ and e₂ = e₁[1 + e₁(m² − 3L)/D₁]; hence M ≈ Ω_t0 and C·Ω_m ≈ Ω_0(h), but neither relation is an identity. The S prefactor limits the leakage of θ-dependence, while the smallness of C·M − Ω_0 and the CMB-fitted track h(θ) control the actual residual. Along the verified 6.1°–30° tracks, C·Ω_m varies by (max−min)/min = 0.11417% (P1) and 0.11924% (P2), whereas Ω_m itself varies by about 32%. This approximate regularity does not survive as observational degeneracy: ℓ_a is empirically almost constant, but R rises by 11.03% and the compressed-CMB χ² rises to 2369.29/2362.06 (P1/P2).

### 2.4. B-Parameters
The B-parameter formalism (Högås & Mörtsell 2021 [Hogas2021]) provides an alternative parametrisation of the bimetric model, where the interaction potential parameters β_n (n = 0,...,4) are re-expressed through observables B_0 and B_1 with direct physical interpretations: B_0 characterises the effective cosmological constant and B_1 the massive graviton contribution. The parameters α and β used throughout this work are reparametrised combinations of the β_n (see Högås & Mörtsell 2021 for the exact mapping). An independent parametrisation of the same model class in terms of (Λ, m_FP, ᾱ), with supernova constraints, was explored by Lüben, Schmidt-May & Weller (2020) [Luben2020b]. The numerical values (Table A1) are computed from the model parameters (θ, m_FP, α, β, Ω_Λ) via the B-parameter definitions in Högås & Mörtsell (2021, Eq. 3.4) [Hogas2021], which express B_0 and B_1 as functions of the interaction potential parameters β_n and the cosmological parameters. B_0 = -40.248, B_1 = 28.230, B_2 = -21.173, B_3 = 21.173, B_4 = -28.206, B_1² = 796.94 (from full-precision value; 28.230² = 796.93 at displayed precision). All B-parameter values are recovered to a precision of 10⁻¹³ in the internal computation, while E²(0) = 1 is recovered to ~10⁻¹⁵, confirming the consistency of the formalism. The round-trip test verifies that the full computational pipeline is self-consistent: starting from the model parameters (θ, m_FP, α, β, Ω_Λ), we solve the quartic for y(z), compute E²(z), reconstruct the quartic coefficients, and recover the original model parameters to a precision of 10⁻¹³. This ensures that no information is lost or corrupted in the forward–backward computation cycle. Note: the displayed B-parameter values are rounded to 3–6 significant figures; exact recovery from these rounded values yields deviations of order 10⁻⁴, not 10⁻¹³. The 10⁻¹³ precision refers to the full-precision internal computation.
### 2.5. Effective Equation of State
The effective dark energy equation of state is w_DE(z) < -1 (phantom) at all redshifts where Ω_DE is significant (i.e., z ≲ 5; at z ≈ 4.9, Ω_DE → 0 and w_DE becomes ill-defined), with w_DE(0) = -1.005 (Table A3). This is a genuine prediction of the bimetric model, arising from the y-dependence of Ω_DE. In the bimetric framework, this phantom regime arises naturally from the massive graviton dynamics and does not require additional phantom fields (Fig. 2).
*Figure 2: Effective dark energy equation of state w_DE*
[^1]: The main text of Högås & Mörtsell (2025) defines y as a_g/a_f, which is inconsistent with their Eq. A10 and with the physical requirement y → 0 at the Big Bang; we follow Eq. A10.

### 2.6. Chi-Squared Calculation and Baseline Parameters
For each probe, χ² = (data - model)^T C^(-1) (data - model), minimised over the relevant free parameters. The resulting per-probe quantities are profile-likelihood diagnostics: they are not integrated marginal posteriors and do not by themselves define evidence or a Bayes factor. For BAO, we profile f = c/(H₀ r_d) by χ² minimisation; for SNe, we profile M analytically (and Ω_m for ΛCDM; the bimetric Ω_m is fixed); for CMB, over h and ω_b (and ω_Λ for ΛCDM; the bimetric Ω_Λ ≈ 0.699 is fixed by the model). The bimetric parameters (θ, m_FP, α, β) are held fixed at the baseline point θ = 6.1°, m_FP/H₀ = 25, α = 0, β = 3, which follows the best-fit configuration of Högås & Mörtsell (2025) for the DESI-era data combination; this choice anchors the comparison to the observationally favoured region of the parameter space. Radiation density: The BAO computation uses Ω_r = 5.38 × 10⁻⁵ (photons only, consistent with the DESI DR2 fiducial cosmology), while the SNe and CMB computations include the neutrino contribution (N_eff = 3.046, Ω_r = 9.08 × 10⁻⁵); this probe-dependent convention follows the fiducial treatment of each dataset, and the effect on Δχ² is < 0.01 for both probes, as Ω_r is negligible at z < 3. Numerical caveat: at the degenerate locus 1+2α+β = 0 (e.g., (α, β) = (−1, 1)) the leading coefficient of the quartic vanishes; for z ≳ 0.76 the solver then falls back to y = 1 (with a stderr warning), reverting E²(z) to its ΛCDM-like form. This point lies off the physical branch (B₁ = 0) and does not affect the canonical results, but grid values there should be read with this caveat. For CMB, absolute χ² and the comparative Δχ² are distinguished explicitly in §5.2; BAO and SNe Δχ² values are relative to their respective ΛCDM best fits. The ΛCDM CMB reference minimises over the three physical parameters (ω_m, h, ω_b), with Ω_r = ω_r/h², Ω_ν = ω_ν/h², Ω_m,total = ω_m/h² and Ω_Λ = 1 − Ω_m,total − Ω_r.
## 3. BAO
### 3.1. Data
We use the DESI DR2 BAO measurements [DESI2025a, DESI2025b], building on the original BAO detection by [Eisenstein2005], consisting of 13 data points spanning redshifts z = 0.295–2.33: D_M/r_d, D_H/r_d, and D_V/r_d. The DESI DR2 covariance matrix (Cobaya release) is block-diagonal, with 2×2 blocks for each tracer pair (D_M, D_H) and a 1×1 block for D_V — there are no cross-tracer correlations. Using exact σ values from the Cobaya matrix yields Δχ² = 0.204. Using rounded σ values from Table IV of the DESI DR2 paper shifts the result by ≈ 10%; the exact value (0.204) is adopted as the baseline. The BAO observables are computed as D_M/r_d = c/H₀ r_d ∫₀ᶻ dz'/E(z'), D_H/r_d = c/H₀ r_d · 1/E(z), D_V/r_d = [z · (D_M/r_d)² · (D_H/r_d)]^(1/3). The parameter f = c/(H₀ r_d) is profiled by χ² minimisation using Brent's method on the interval [20, 80], which spans the physically expected range f ≈ 25–35 (for r_d ≈ 147 Mpc and H₀ ≈ 68–74 km/s/Mpc) with generous margins.
### 3.2. Results
The bimetric model yields χ²_BAO = 10.472 vs. ΛCDM χ²_ΛCDM = 10.284, giving Δχ² = 0.188. The formal one-parameter tail value p = 0.664 is retained only as a descriptive reference, not as a model-selection significance. Wilks-regularity conditions are not automatic here; in particular, at θ = 0 the parameter m_FP becomes unidentifiable. No detection or preference probability is inferred from this tail value. In the reported baseline scans, BAO is insensitive to the four varied bimetric parameters at current data precision, with partial degeneracy under profiling of the distance ratio f:
- (α, β) scan (20 combinations): Δχ² ≈ 0.188 — flat across the entire grid.
- m_FP/H₀ scan (10 ≤ m_FP/H₀ ≤ 10⁵, four orders of magnitude): Δχ² varies only weakly across the grid (remaining far below the 3.84 threshold), with the baseline value 0.188 at m_FP/H₀ = 25.
- θ scan (block-diagonal covariance): the self-consistent scan (Ω_m refit at each θ) is flat — Δχ² decreases mildly from 0.204 at θ = 6.1° to 0.198 at θ = 31.5°: the profiling of f absorbs the θ-dependence even at large angles. Only a fixed-Ω_m variant (no refit) shows growth, reaching Δχ² = 56.1 at θ = 31.5° — this variant overstates the sensitivity, since Ω_m in the model is set by the quartic constraint, not held fixed. Our θ = 6.1° is deep in the blind regime. Under the ΛCDM null (10⁴ synthetic DESI DR2 realisations drawn from the full covariance, with f refit for both models), Δχ² ≥ 0.188 occurs with probability p = 0.689: the BAO penalty is statistically unremarkable.
The Δχ² decomposition reveals that the dominant contribution comes from the Ω_m shift (+0.682), which is partially compensated by the bimetric E²(z) correction (-0.494), yielding the net Δχ² = 0.188. The pure bimetric correction at fixed ΛCDM Ω_m contributes only Δχ² = 0.021; this is a standalone diagnostic and is not an additive term in the decomposition (0.682 - 0.494 = 0.188, while 0.682 - 0.494 + 0.021 = 0.209 ≠ 0.188). The bimetric E²(z) correction compensates for 72% of the Ω_m penalty (-0.494 out of +0.682). The 70× ratio of SNe pure correction (1.49) to BAO pure correction (0.021) arises because BAO observables are dimensionless ratios D_M/r_d and D_H/r_d, where the profiling of f = c/(H₀ r_d) absorbs the approximately constant component of the bimetric E²(z) correction, leaving only a residual z-dependent signal of order 0.021. SNe observables (luminosity distances D_L) have no analogous scale-profiling degeneracy: the absolute scale is set by H₀ through M, but the shape of D_L(z) directly reflects E²(z), so the full 1% correction contributes 1.49.
## 4. SNe Ia
### 4.1. Data
We use the Pantheon+ sample (Brout et al. 2022 [Brout2022]; Scolnic et al. 2022 [Scolnic2022]), consisting of 1701 light curves from 1550 unique supernovae. After applying the standard z > 0.01 cut (removing 113 low-redshift calibrators), 1588 light curves remain in the redshift range z = 0.01–2.26. Högås & Mörtsell (2025) employ a stricter z > 0.023 cut. The choice of z > 0.01 shifts the best-fit ΛCDM Ω_m to higher values, potentially enhancing the apparent tension with the bimetric prediction Ω_m = 0.304. We assess the robustness of the SNe sensitivity to this choice via the σ_int scan (§4.2) and the full-covariance analysis (§4.2). The distance modulus is μ = 5log₁₀(D_L) + 25, where D_L = (1+z) D_M, and M is treated as a free nuisance parameter and profiled jointly with Ω_m for the ΛCDM fit; for the bimetric model, Ω_m is predicted by the quartic constraint (§2.3) and held fixed, and only M is profiled.
Covariance treatment. We report two statistically consistent estimates. (i) Diagonal approximation (m_b_corr_err_DIAG): neglects off-diagonal systematic correlations. This approximation yields Ω_m^best = 0.374 for ΛCDM, a 2.2σ shift from the Pantheon+ full-covariance result (Ω_m = 0.334 ± 0.018 [Brout2022]), attributable to the neglected correlations between low-z and high-z supernovae. This shift directly affects Δχ²_SNe. (ii) Full covariance (STAT+SYS, 1588 × 1588 after the cut, with analytical profiling of M): the ΛCDM best-fit moves to Ω_m = 0.347 and Δχ²_SNe decreases from 15.57 to 6.50 (§4.2). The result is robust to the low-redshift cut: with the Högås & Mörtsell cut z > 0.023 (1367 light curves), the full-covariance values are Δχ²_SNe = 5.79 (pure 0.88) and Ω_m^best = 0.346. The diagonal value 15.57 should therefore be interpreted as an upper bound on the SNe tension. A σ_int scan with the diagonal covariance (Table A6) provides an intermediate, conservative estimate (11.58 at σ_int = 0.12 mag).
### 4.2. Results
Diagonal covariance. The bimetric model yields Δχ²_SNe = 15.57 for the baseline parameters, using the standard comparison: ΛCDM best-fit (Ω_m = 0.374, M profiled) vs. bimetric (Ω_m = 0.304 from quartic, M profiled), χ²/dof = 0.455. SNe Ia are not blind to the bimetric model:
- Standard comparison (diagonal): Δχ² = 15.57 (Ω_m^LCDM = 0.374 vs. Ω_m^bim = 0.304).
- Pure bimetric (same Ω_m = 0.304): Δχ² = 1.49 — the E²(z) correction alone is sub-detectable.
- (α, β) scan: Δχ² = 13.9–15.6. m_FP/H₀ scan (≥ 5): Δχ² = 15.57. θ scan: Δχ² = 13.9–15.6 (flat, but non-zero).
Full covariance (STAT+SYS). Repeating the analysis with the full 1588 × 1588 Pantheon+ STAT+SYS covariance (analytical profiling of M) gives: ΛCDM best-fit Ω_m = 0.347 (χ²_min = 1433.06, χ²/dof = 0.903); bimetric χ² = 1439.56; Δχ² = 6.50 (standard), and Δχ²_pure = 0.96. The full covariance thus reduces Δχ² by ~58%: the systematic correlations pull the ΛCDM best-fit Ω_m from 0.374 to 0.347 — closer to the bimetric prediction (0.304) and to the BAO/CMB-preferred range — which shrinks the Ω_m mismatch that dominates Δχ². The decomposition mirrors the diagonal case: of the 6.50, 5.54 is attributable to the Ω_m shift and 0.96 to the pure E²(z) correction. The qualitative conclusion (SNe sensitive through the Ω_m prediction) is unchanged. Square-root-Δχ² values are not calibrated as standalone σ preferences and are not reported as such.
The Δχ² decomposition (diagonal case) reveals: Δχ² = -2 × cross-term + model-diff χ², where the cross-term measures the correlation of data residuals with the model difference, and model-diff χ² measures the centred shape difference. At Ω_m^LCDM = 0.374: model-diff = 15.64, cross-term = 0.03, giving Δχ² = 15.57. The Ω_m shift dominates: the bimetric E²(z) correction is negative (≲ 1%), which worsens the fit when Ω_m^bim < Ω_m^(LCDM best). For BAO, the same correction compensates (Ω_m^bim = 0.304 > Ω_m^(BAO best) = 0.297), but for SNe the direction is reversed. With the full covariance, this picture persists in weakened form: model shape differences contribute only 0.96, while the remaining 5.54 reflects the compressed Ω_m offset (0.347 - 0.304 = 0.043).
Intrinsic dispersion. To assess the robustness of the SNe sensitivity to the treatment of intrinsic scatter, we additionally compute Δχ²_SNe including a standard intrinsic dispersion σ_int = 0.12 mag, added in quadrature to the diagonal covariance (σ²_total = σ²_diag + σ²_int). This reduces Δχ²_SNe from 15.57 to 11.58 (Δχ²_pure from 1.49 to 1.11), while Ω_m^best remains 0.374. Even at σ_int = 0.20 mag (an upper bound), Δχ²_SNe = 8.06 — SNe Ia remain a sensitive discriminator. The qualitative conclusion is unchanged. We do not add σ_int on top of the full STAT+SYS covariance: the intrinsic dispersion is already included in its diagonal (standard Pantheon+ construction), and adding it again would double-count the same physical effect. The σ_int scan (diagonal) is tabulated in Table A6.
The Ω_m = 0.374 obtained from the diagonal covariance is shifted 2.2σ from the Pantheon+ full-covariance value (0.334 ± 0.018 [Brout2022]). With full covariance, Ω_m^(LCDM best) = 0.347, closer to the bimetric prediction, reducing Δχ²; the Ω_m shift decreases from 0.070 to 0.043. This is a quantitative, not qualitative, change: SNe remain more sensitive to the bimetric Ω_m prediction than BAO under the stated reference. The fixed-benchmark Δχ² = 6.50 is not assigned a universal 95% model-selection threshold.
## 5. CMB
### 5.1. Data: Compressed CMB Likelihood
We use the compressed Planck 2018 + ACT DR6 likelihood from Bansal & Huterer (2025) [Bansal2025], which provides a 3-parameter data vector: R = 100·√(Ω_b h² + Ω_cdm h² + Ω_ν,m h²) · D_M(z_*)/c (shift parameter, following Bansal & Huterer 2025, Eq. 4). Tensor modes in ghost-free bimetric gravity propagate at the speed of light (c_gw = c), so the model trivially satisfies the GW170817 multimessenger bound. Massive neutrinos (Σm_ν = 0.06 eV) are evolved as nonrelativistic matter at z ≲ 130 and as radiation above, following the convention of the compressed likelihood (Bansal & Huterer 2025, Appendix C) — convention P1 in the terminology of §7.8 (item 5); the R prefactor follows their Eq. 4 exactly, ℓ_a = π D_M(z_*)/r_s(z_*), Ω_b h² = 0.022371; the observed data vector is (R, ℓ_a, ω_b) = (1.7504, 301.77, 0.022371), with the 3×3 covariance matrix (in units of 10⁻⁸): C = [1559.83, −1325.41, −36.45; −1325.41, 714691.80, 269.77; −36.45, 269.77, 2.10] × 10⁻⁸. The sound horizon is r_s(z_*) = ∫_{z_*}^∞ c_s(z)/H(z) dz, c_s(z) = c/√(3(1+R_b(z))), R_b(z) = 3Ω_b / [4Ω_γ(1+z)], where Ω_γ h² = 2.469 × 10⁻⁵ and z_* = 1090. The neutrino effect on r_s is primarily through H(z), where the total radiation density Ω_r includes the neutrino contribution. The sound speed c_s = c/√(3(1+R_b)) depends only on the photon and baryon densities and does not require a neutrino correction (Hu Sugiyama 1995 [HuSugiyama1995]). The parameter Ω_Λ ≈ 0.699 is held fixed at the H&M best-fit anchor (not varied in the fit); h and ω_b are the two free parameters in the CMB fit.
Applicability to modified gravity: The compressed CMB likelihood is derived under GR and the Planck Collaboration (2015, Results XIV) [Planck2015XIV] warns that it may be inapplicable to modified gravity models where perturbation growth differs from GR. We use the likelihood only conditionally. Vainshtein screening suppresses the background modification at z_* = 1090, where H(z_*) ≫ m_FP, but screening alone does not validate the transfer of a GR-calibrated likelihood to bimetric perturbations. Lüben et al. (2020) [Luben2020] used CMB data in their bimetric analysis, justifying it by the same Vainshtein screening argument. We stress that this screening argument predates the scalar-sector analysis of §6: the gradient instability at z > z* (X(1090) < 0) is a background-driven property of the finite branch that Vainshtein screening does not remove. The perturbation-level applicability of the compressed likelihood is therefore conditional (§6.4), and our CMB results should be read as a background-level geometric diagnostic with this caveat. **Warning — perturbations not included.** The compressed likelihood probes the background geometry only. Beyond the shift parameters, bimetric gravity modifies the growth of structure and the decay of gravitational potentials throughout the unscreened region (z ≲ 12), sourcing an Integrated Sachs–Wolfe signal accumulated along the line of sight; a full Boltzmann treatment (CAMB/CLASS with the bimetric stress-energy tensor) could therefore constrain — or exclude — the model far more strongly than the background χ² suggests, and is the single most important follow-up. CMB E²(z): The full bimetric E²(z) is used for the D_M(z_*) integral. The ΛCDM E²(z) approximation, while adequate at θ = 6.1° (correction < 0.3%), underestimates the absolute CMB χ² by up to 3.6× at intermediate θ and is not used.
### 5.2. Results
The comparative CMB reference is a flat ΛCDM model with three free physical parameters (ω_m, h, ω_b), rather than the former two-parameter fit at fixed matter and dark-energy anchors. We use Ω_r = ω_r/h², Ω_ν = ω_ν/h², Ω_m,total = ω_m/h², Ω_m,nonrel = Ω_m,total − Ω_ν and Ω_Λ = 1 − Ω_m,total − Ω_r. Three compressed data points and three fitted parameters give χ²_min ≈ 0 and zero residual DOF; therefore no reference goodness-of-fit p-value is defined. The independently reproduced references are:

```
Convention   χ²_ref      h          ω_m         Ω_m,total   Ω_m,nonrel   Ω_Λ
P1           < 10⁻¹²     0.711386   0.159862    0.315889     0.314616      0.684029
P2           < 10⁻¹²     0.672159   0.142957    0.316419     0.314993      0.683489
```

Throughout this paper Δχ²_CMB = χ²_CMB,bim − χ²_CMB,ΛCDM, with the reference computed in the same neutrino convention. Since χ²_ref ≈ 0, the comparative and absolute values coincide to the displayed precision; the former negative baseline offsets are superseded. The bimetric baseline has two fitted parameters (h, ω_b), hence DOF = 1. It gives χ² = 0.0409707 at h* = 0.724201 in P1 (p = 0.840) and χ² = 0.0771240 at h* = 0.684662 in P2 (p = 0.781). Both baselines are consistent with the compressed data, but neither the absolute p-values nor the small positive Δχ² constitutes a model-selection significance. The discriminatory power comes from the θ-dependence:

- Baseline (θ = 6.1°): Δχ² = 0.04097 (P1) and 0.07712 (P2), with χ²_ref ≈ 0 in each convention.
- (α, β) scan (20 combinations): absolute χ² = 0.014–2.14, hence Δχ² = 0.014–2.14 against the free P1 reference. The predicted Ω_m is independent of (α, β) across the entire grid (0.3029 at all 20 points — convention (ii) of the Appendix A note; |∂Ω_m/∂α| ≈ 4 × 10⁻⁹, |∂Ω_m/∂β| ≈ 10⁻¹², Appendix B.2), because the de Sitter attractor pins y(0) to 1 − 4.9 × 10⁻⁴ and the (α, β)-dependent terms in Ω_DE(0) enter only at O((1−y(0))²). Consequently, the SNe Ω_m tension is not alleviated across this scanned grid; this is a restricted numerical and local-analytic result, not a global proof over all (α, β) (§7.5). The m_FP/H₀ scan (5 ≤ m_FP/H₀ ≤ 100) gives absolute and comparative χ² = 0.001–0.36 in the same reference convention.
- θ scan: in the unified independently reproduced vintage, P1 gives χ² = 0.04097 at 6.1° and 2369.29001 at 30°; P2 gives 0.07712 and 2362.06129. Thus Δχ² = 0.04097 → 2369.29001 (P1) and 0.07712 → 2362.06129 (P2). The full tables are given in Appendix A5. At θ = 30° with m_FP/H₀ = 25, the static Higuchi bound m_FP² > 2Ω_Λ is satisfied (625 > 1.4), although the angle is strongly penalised within the compressed-CMB geometric diagnostic.

The split test isolates the source of the θ penalty in the D_M integral. The residual E²_bim − E²_screened is set to zero only at z > 12 (`zero_high`) or only at z < 12 (`zero_low`), while the screened sound horizon is retained in every mode. In P1, for example:

```
θ       full χ²     zero_high χ²   zero_low χ²
10°     9.35910     7.58224        0.77763
30°     2369.29001  2140.90011     73.64757
```

At 30° the mode retaining the late-time residual therefore retains 90.36% of the full χ², whereas the mode retaining only the early-time residual retains 3.11%. This is a diagnostic of the distance integral, not a full CMB-spectrum calculation and not evidence for an ISW or physical-θ mechanism.

The CMB-only sensitivity to θ, after isolating the BAO and SNe contributions, is presented for the first time in this work (Fig. 3). Högås & Mörtsell (2025) obtained θ constraints from the combined fit but did not decompose the constraining power by probe. The θ grid (6.1°–30°) spans the posterior support of Högås & Mörtsell (2025) and the angles that would reproduce the externally preferred Ω_m values; the m_FP/H₀ scan spans 5–100 (≈1.3 decades).

### 5.3. Summary of χ² and Δχ² Scan Results
Table 2: Summary of scan results. BAO and SNe entries are Δχ²; both absolute and comparative CMB values are shown. The free three-parameter ΛCDM reference has χ²_ref ≈ 0, so the last two CMB columns coincide to the displayed precision.
```
Scan              BAO Δχ²      SNe Δχ²                CMB χ² (P1)   CMB Δχ² (P1)
(alpha,beta), 20  0.188        13.9--15.6 (diag)      0.014--2.14    0.014--2.14
combos
m_FP/H0           0.188 (weak    15.57 (diag)           0.001--0.36    0.001--0.36
                  dependence)
theta (baseline   0.204 (bd)   15.57 (diag); 6.50     0.04097        0.04097
6.1deg)                        (full cov)
theta (up to      flat ->0.198  13.9--15.6 (diag)       -> 2369.29     -> 2369.29
31.5deg)          (56.1 fixed Om)
Verdict           Blind        Sensitive               Sensitive to theta
```
P2 gives 0.07712 → 2362.06129 over the same θ grid, with the free P2 reference. P1 and P2 references are not interchangeable.
![Figure 3](figures/Fig3_theta_scan.png)

*Figure 3: Per-probe Δχ² as a function of the mixing angle θ. BAO (blue) stays flat (Δχ² ≈ 0.20, self-consistent block-diagonal scan; only the fixed-Ω_m variant grows, reaching 56.1 at θ = 31.5°), while SNe Ia (red) are offset at Δχ² ≈ 15 (diagonal; 6.5 with full covariance). The CMB curve (green) uses the free three-parameter P1 ΛCDM reference (χ²_ref ≈ 0), rising from 0.04097 at baseline to 2369.29001 at θ = 30°. A symmetric-logarithmic vertical scale retains the small positive baseline; the dotted line shows Δχ² = 3.84 (95% reference for 1 DOF where applicable).*
## 6. The Scalar Sector: X(z), z*, and the Growth Factor Γ

**Convention line (fixed for this paper):** perturbations Φ ~ exp(iωN) with N = ln a (Akrami et al. 2015, Eq. 37 [Akrami2015]); instability when X < 0, i.e. ω² < 0; the growth integral is evaluated in the canonical 1/E convention (the dimensionally inconsistent 1/E² variant is not used anywhere in this work). Numerical values below come from the project pipeline and were reproduced by an external blind replica (GLM, agreement to 8 decimal places on the reproducible code value of z*, see Appendix C).

### 6.1. X(z) and the instability boundary
Background-level viability is necessary but not sufficient: a cosmological model must also be free of ghost and gradient instabilities. We evaluate the squared eigenfrequency of the scalar perturbation kinetic sector, X(z) ≡ (aH/k)²ω², in the sub-horizon WKB limit, using Eq. 38 of Akrami et al. (2015) [Akrami2015] (the eigenfrequencies were originally presented in Könnig et al. 2014 [Konnig2014]) with the mapping between the Högås–Mörtsell and Akrami conventions derived from the rescaling symmetry of the bimetric action (their Eq. 43): in the c = 1 gauge, α_A = tan θ, and the consistency condition m²β_n^(A) = B_n holds identically for n = 1, 2, 3 with m² = m_FP² sin²θ = 7.0575, giving (β_1, β_2, β_3)^(A) = (4, −3, 3) — the pure Weinberg parameters. The mapping identity was verified numerically to 10⁻¹³.

The function X(z) governs the scalar sector: X < 0 signals the gradient instability, with the growing mode exp(iωN) exponentially amplified while the unstable window stays open. At the baseline parameters (θ = 6.1°, m_FP/H₀ = 25, α = 0, β = 3) the canonical evaluation gives:

- X(0) = +0.999996 (full-precision canonical value 0.9999955766);
- X(5) = +0.839993; X(10) = −0.465722;
- X(1090) = −1.163845.

The sign structure is: X > 0 throughout the observable late universe, crossing to X < 0 at the boundary z* = 7.880536 (the superseded office-pipeline value 7.8836 is documented in Appendix C). The epoch z* lies in the reionisation era — well before recombination (z = 1090). The growing regime is therefore confined to the epoch between z* and recombination: the window opens at z* and extends back into the past toward z = 1090; it does not extend to the present epoch, where X > 0 and the linear sector is perturbatively stable.

This localisation is quantitatively consistent with the finite-branch behaviour of Akrami et al. (2015) [Akrami2015] and with the no-go analysis of Könnig (2015) [Konnig2015]: the gradient instability sets in at early times on the finite branch; at the baseline point the crossing lies at z* = 7.880536, where H ≈ 14.6 H₀ — below m_FP = 25 H₀, i.e. the instability boundary is distinct from the Vainshtein screening boundary (z ≈ 12, §7.1). Our result locates the boundary quantitatively for the full (all B_n ≠ 0) model at the baseline point.

**Three distinct redshifts.** The boundary z* = 7.880536 is a perturbative root: it is defined by X(z) = 0 with X(z) ≡ (aH/k)²ω² (Eq. (38) of [Akrami2015], convention Φ ∝ exp(iωN)), and is not a background-level crossing. It must not be confused with (i) the background zero z₀ of Ω_DE(z) — z₀ = 4.882 in the radiation-only convention of this section and z₀ = 4.890 in the neutrino-closure convention of Appendix A — at which the background Higuchi sensor m²_eff = 2Ω_DE changes sign, nor with (ii) the Vainshtein screening boundary z_V = 11.73 (≈ 12), defined by H(z_V) = m_FP at m_FP/H₀ = 25 (§7.1). At the baseline point the three are ordered z₀ < z* < z_V and are conceptually independent: the first two characterise the background branch, the last the screening regime, while z* alone governs the linear scalar sector.

**Caveat.** Strictly, linear perturbation theory is not formally applicable during the unstable epoch z > z*: the linearised equations lose their domain of validity once a gradient instability switches on. The statements above should be read as locating the instability boundary, not as evolving perturbations through it. Nonlinear evolution of a homogeneous spherical overdensity in this theory remains stable (Högås, Torsello & Mörtsell 2020 [Hogas2019stab]), so the instability may be bounded by nonlinearities — deciding this requires a full Boltzmann-level treatment, which does not currently exist for Hassan–Rosen bimetric gravity (the closest existing computations target low-redshift observables — fσ₈, skewness S₃, ISW cross-correlations — without a full Boltzmann treatment, e.g., Bassi et al. 2023 [Bassi2023]).

### 6.2. The growth integral
Growth of the unstable mode is quantified by the accumulated exponent

ln Γ(k) = (k / H₀) · I₁,   I₁ = ∫_{z*}^{1090} √(|X(z)|) / E(z) dz,

with the canonical 1/E convention (dimensionally consistent: the factors (1+z) from dN = −dz/(1+z) and aH = H/(1+z) cancel exactly, leaving k dz/(H₀E)). The integrand carries the modulus √|X|/E, since X < 0 throughout the growing window. The integration runs from z* back into the past, to z = 1090 — the last-scattering boundary of the integration domain (X remains negative at 1090; the integral is truncated at the epoch where the compressed CMB likelihood applies). The window does not extend to the present epoch.

The canonical value of the integral is I₁ = 0.945210110 (evaluated with the reproducible code value z* = 7.880536; verified by an independent canonical-code run; the full Γ-table below is pure arithmetic on I₁ with H₀(h) = 100h/299792.458 Mpc⁻¹ and is not recomputed per row). The two values of h are the reference values of the analysis: h = 0.7 (baseline) and h = 0.678 (the neutrino-closure reference value of the Appendix A note):

```
k [Mpc^-1]   h      k/H0      ln Gamma   Gamma          log10 Gamma
0.001        0.700  4.2827    4.0481     57.29          1.76
0.002        0.700  8.5655    8.0962     3.28 x 10^3    3.52
0.003        0.700  12.8482   12.1443    1.88 x 10^5    5.27
0.01         0.700  42.8275   40.4810    3.81 x 10^17   17.58
0.05         0.700  214.1375  202.4049   8.00 x 10^87   87.90
0.10         0.700  428.2749  404.8098   6.41 x 10^175  175.81
0.001        0.678  4.4217    4.1795     65.33          1.82
0.002        0.678  8.8434    8.3589     4.27 x 10^3    3.63
0.003        0.678  13.2652   12.5384    2.79 x 10^5    5.45
0.01         0.678  44.2172   41.7945    1.42 x 10^18   18.15
0.05         0.678  221.0859  208.9726   5.70 x 10^90   90.76
0.10         0.678  442.1718  417.9452   3.25 x 10^181  181.51
```

Consistency check: k = 0.05 Mpc⁻¹, h = 0.7 → log₁₀Γ = 87.90. Derived values: Γ(0.001 Mpc⁻¹) ≈ 57, Γ(0.003 Mpc⁻¹) ≈ 2 × 10⁵.

### 6.3. Baseline linear-sector diagnostic
At the baseline point, growth of 10^87.9 over the unstable window means that the adopted linear/WKB scalar-sector description is unviable: a seed perturbation in that formal window is amplified by ~88 orders of magnitude, and linear perturbation theory fails long before structure formation can be described. This is not a failure of the entire bimetric theory or of every parameter point — it is the finite-branch gradient instability behaving as the literature predicts (Akrami et al. 2015 [Akrami2015]; Könnig 2015 [Konnig2015]), now with a precise number attached.

The reviewer's fork — whether the instability is softened on the largest scales — has a quantitative answer: on near-horizon / largest scales, k ≈ 0.0011–0.0023 Mpc⁻¹, Γ falls to the 10^2–10^4 range. The separation between the baseline point and the survival window is one of 84–86 orders of magnitude — of degree, not of kind: even in the most favourable corner of the fork, the accumulated growth is enormous, though orders of magnitude below the sub-horizon value.

For contrast, the dynamical Higuchi bound — the generalisation of the static bound m_FP² > 2Ω_Λ (Högås & Mörtsell 2021, Eq. 5.11) evaluated along the expansion history — is satisfied over the entire history at the baseline parameters, with a minimum margin of 892×: the theory is not ruled out by ghost excitation. The pathology is specifically the finite-branch scalar gradient instability.

### 6.4. Required viability tests
Within the adopted linear/WKB diagnostic, the baseline scalar sector is unviable (Γ = 10^87.9 at k = 0.05 Mpc⁻¹). Possible ways to obtain a viable regime include parameter regions where z* is pushed outside the observable window or a nonlinear saturation mechanism; neither is established by this paper. All structure-formation statements — late-ISW, S₈, and growth — therefore remain conditional background-level statements. A certified viability map and an accepted joint analysis remain necessary; no accepted production MCMC result is claimed in this draft.

## 7. Discussion
### 7.1. Vainshtein Screening: Per-Probe Mechanism
The blindness of BAO and the sensitivity of SNe and CMB arise from the interplay of two distinct mechanisms: the de Sitter attractor (y → 1) that suppresses the bimetric correction at late times, and the Vainshtein screening (Lüben et al. 2020 [Luben2020]) that operates in the early universe. The separation of these mechanisms across redshift ranges supports the proposed interpretation within the tested baseline and grids; it does not prove the same pattern over the full model class.
The critical Hubble rate is H_c = m_FP (Lüben et al. 2020). For orientation, the Vainshtein radii implied by m_FP/H₀ = 25 are r_V(☉) ≈ 11 pc and r_V(Milky Way) ≈ 100 kpc; the GW170817 bound on the gravitational-wave speed is ~9 orders of magnitude weaker than would be needed to exclude this class, so screening neither rescues nor undermines the cosmological conclusions — it delineates their regime of applicability. For H > m_FP (early universe, z ≳ 12 at m_FP/H₀ = 25), Vainshtein screening is active and GR is restored. For H < m_FP (late universe, z ≲ 12), screening is inactive and the bimetric correction is present. In this regime, GR is approached not through Vainshtein screening but through the de Sitter attractor y → 1 (Lüben et al. 2020), which makes the correction small but non-zero.
BAO measures D_M/r_d and D_H/r_d at z < 2.33 — the late universe where Vainshtein screening is inactive. The bimetric correction to E²(z) is present but intrinsically small (≲ 1.1% within the BAO range), because the dynamical evolution drives y → 1 at late times as a de Sitter attractor (y(0) = 0.99951), making Ω_DE ≈ Ω_Λ with only a tiny correction. The profiling of f = c/(H₀ r_d) further absorbs residual signal, because BAO observables are dimensionless ratios. The combination renders Δχ² ≈ 0.
SNe Ia measure D_L(z) = (1+z) D_M(z) at z < 2.26 — the same late-time regime where screening is inactive. The H(z) integral accumulates only the small unscreened correction (≲ 1.2%), and the pure bimetric correction at fixed Ω_m gives Δχ² = 1.49 (diagonal) / 0.96 (full covariance) — sub-detectable. However, SNe Ia are sensitive to the bimetric model through the Ω_m shift: the quartic constraint predicts Ω_m = 0.304, while ΛCDM best-fit to Pantheon+ gives Ω_m = 0.374 (diagonal) / 0.347 (full). The resulting Δχ² = 15.57 / 6.50 is dominated by this Ω_m mismatch, not by the E²(z) correction itself. Unlike BAO, there is no analogous r_d-scale profiling degeneracy: D_L depends directly on H₀ and Ω_m through the luminosity distance integral. The negative bimetric E²(z) correction shifts the model in the wrong direction relative to the SNe-preferred Ω_m, amplifying rather than compensating the Ω_m discrepancy. This contrasts with BAO, where Ω_m^bim = 0.304 > Ω_m^(BAO best) = 0.297 and the correction partially compensates.
CMB measures D_M(z_*) — the comoving distance integrated from z = 0 to z_* = 1090 — and the sound horizon r_s(z_*) = ∫_{z_*}^∞ c_s(z')/H(z') dz' (from z_* to ∞). The integral D_M(z_*) traverses two regimes: the unscreened late universe (z ≲ 12, where H < m_FP and the bimetric correction is present) and the screened early universe (z ≳ 12, where H > m_FP and the correction is suppressed). Although the correction in the unscreened regime is small (≲ 1.2%), the unscreened region (z ≲ 12) contributes 72% of D_M(z_*) (Table A4; computed as ∫₀¹² dz/E(z) / ∫₀¹⁰⁹⁰ dz/E(z) ≈ 0.72 for the bimetric model with baseline parameters), so the correction accumulates and produces a detectable shift in D_M(z_*), which the compressed CMB likelihood registers as a change in absolute χ² and, away from the baseline, as a large positive Δχ² relative to ΛCDM. The sound horizon r_s(z_*), being an integral from z_* to ∞, lies entirely in the screened regime and is less affected. This asymmetry between the screened early universe and the unscreened late universe is the fundamental reason CMB is a sensitive discriminator among the three probes, alongside SNe which are sensitive through the Ω_m prediction (Table 3).
Table 3: Per-probe blindness and sensitivity mechanism summary.
```
Probe    Bimetric        Nuisance treatment  Source of sensitivity/blindness
         correction to                          
         E2(z)                                 
BAO      ~< 1.1% (z<=2.33) Yes: f = c/(H0 rd)   Blind: small correction +
                                             f-profiling absorbs
                                             residual
SNe Ia   ~< 1.1% (z<=2.26) Yes: M profiling    Sensitive: Omega_m mismatch
                                             (0.304 vs 0.374 diag / 0.347
                                             full) dominates (~90% diag,
                                             ~85% full)
CMB (D_M) ~< 1.2% (z ~< 12, No               Sensitive: 72% of D_M(z*) from
          unscreened)                          unscreened region accumulates
                                             correction
CMB (r_s) Suppressed (z >=  No                Blind: r_s integral entirely in
          1090, screened)                      screened regime
```
> Earlier exploratory joint-posterior diagnostics used older conventions and are not accepted P3 posteriors. They are omitted from the evidentiary basis of this public draft; no joint parameter constraint is carried into the developing exact-Fermi–Dirac analysis.

### 7.2. Phantom Dark Energy Phase and the CPL Comparison with DESI
The effective equation of state w_DE(z) < -1 (phantom) at all redshifts where Ω_DE is significant (z ≲ 5) is a genuine prediction of the bimetric model. BAO is insensitive to this phase because the bimetric correction to E²(z) remains small (≲ 1.1% within the BAO redshift range (z ≤ 2.33), rising to ~1.08% at z ≈ 3 and ~1.11% at z ≈ 5 just beyond the BAO range, for context). SNe are sensitive to the bimetric model not through the E²(z) correction directly (Δχ²_pure = 1.49 / 0.96), but through the Ω_m prediction (Δχ² = 15.57 / 6.50). The CPL approximation (w_0 ≈ -0.861, w_a ≈ -0.500) fitted on z [0, 3] does not capture the phantom nature (w_DE(0) = -1.005); Ω_DE(z) crosses zero at z ≈ 4.9 and tends to the negative asymptote Ω_Λ − sin²θ·m_FP²(1 + α + β/3) = −13.42 (at baseline α = 0) (Ω_DE(10) ≈ −3.05), making CPL inapplicable beyond z ≈ 4. The discrepancy arises because the effective w_DE(z) in the bimetric model has a complex redshift dependence — with Ω_DE(z) changing sign at z ≈ 4.9 and becoming negative at higher z — that a 2-parameter CPL fit cannot capture.

For context, Fig. 4 places our CPL fits of the bimetric w_DE(z) alongside the official DESI DR2 CPL results ([DESI2025a], Table V): DESI+CMB (w_0 = -0.42 ± 0.21, w_a = -1.75 ± 0.58), DESI+CMB+Pantheon+ (-0.838 ± 0.055, -0.62 +0.22/-0.19), DESI+CMB+Union3 (-0.667 ± 0.088, -1.09 +0.31/-0.27), and DESI+CMB+DESY5 (-0.752 ± 0.057, -0.86 +0.23/-0.20). Our fits — (w_0, w_a) ≈ (-0.950, -0.252) on the archived 60-row table over z ≤ 2.3 and ≈ (-0.942, -0.267) on a dense grid; independently reproduced by three external codes (Alice AI: -0.9423, -0.2671; DeepSeek: -0.9417, -0.2678; Qwen: -0.9418, -0.2677) — sit visually close to the ΛCDM point (-1, 0); no σ distance is assigned because the full DESI covariance is not used. Two cautions apply. First, the error crosses in Fig. 4 are 1σ without the w_0–w_a covariance, which is substantial for DESI; proper contours require the published Fisher information. Second, the quadrant occupied by a CPL fit of our w_DE(z) depends on the fit range (e.g., z ≤ 0.5 yields w_a ≈ -0.035): the quadrant is an artefact of the approximation, not physics. The robust characteristic is the combination w_0 + w_a ≈ -1.21, common to the fits on the BAO-range tables (the archived z ≤ 2.3 table and the dense grid) and to their external reproductions.

*Figure 4: CPL (w₀, w_a) plane — bimetric fits vs official DESI DR2 combinations.*
### 7.3. Combined Δχ² and the Dominance of SNe
If the three probes are combined, the total Δχ² is dominated by SNe. To keep the combination internally consistent, we use the free three-parameter CMB reference of §5.2 and do not mix covariance choices across probes:
```
Scheme                    CMB     BAO      SNe      CMB Δχ²   Sum Δχ²
All-diagonal              P1      0.204    15.57    0.04097   15.82
All-full                  P1      0.188    6.50     0.04097   6.73
σ_int = 0.12, diagonal    P1      0.204    11.58    0.04097   11.83
All-diagonal              P2      0.204    15.57    0.07712   15.86
All-full                  P2      0.188    6.50     0.07712   6.76
```
The unrounded P1 sums are 15.8191291, 6.7255449 and 11.8258438 respectively. "All-diagonal" uses the block-diagonal BAO result and diagonal SNe covariance; "all-full" uses the full 13×13 BAO covariance and full Pantheon+ STAT+SYS covariance. A hybrid full-BAO/diagonal-SNe combination is not quoted. These sums are descriptive arithmetic under the stated references, not a certified joint likelihood or model-selection statistic; no p-values are assigned. Shared nuisance parameters would require an explicit joint treatment.

SNe supply the dominant positive baseline contribution through the Ω_m mismatch, BAO contributes little, and the new CMB baseline term is small and positive rather than a negative offset. At larger θ the CMB term grows rapidly. In the combined analysis of Högås & Mörtsell (2025), both SNe and CMB drive the constraints on θ, while BAO contributes little.
### 7.4. Implications for Future Surveys
The bimetric correction to E²(z) at BAO redshifts is small (≲ 1.1%) but non-zero — Vainshtein screening is inactive at z < 12, so the correction is present, not suppressed. The current BAO blindness therefore reflects the limited precision of existing data (~1–2%) combined with the profiling of f, not a fundamental theoretical limitation. Future sub-percent BAO surveys (DESI DR3, Euclid) could in principle detect the z-dependent component of the 1% correction, though the f-profiling absorbs any approximately constant (in z) component, leaving only the residual z-dependent signal as potentially detectable. Probes that integrate over longer paths, such as CMB, or which directly sample the transition region — 21 cm cosmology, including the global 21-cm signal at z ≈ 15–30 — remain the most promising avenue for constraining θ.
### 7.5. Ω_m and H₀
We note that the predicted Ω_m = 0.304 is in mild tension with independent measurements: Planck 2018 gives Ω_m = 0.315, and the Pantheon+ full-covariance best-fit gives Ω_m = 0.334 ± 0.018. This tension is itself a constraint on the bimetric model at the baseline parameter point. Two restricted results sharpen the tension at the baseline and along the tested tracks; they do not constitute a global no-go statement. First, the (α, β) scans (§5.2, Table A1b) show that the predicted Ω_m is independent of the potential parameters to |∂Ω_m/∂α| ≈ 4 × 10⁻⁹, |∂Ω_m/∂β| ≈ 10⁻¹²: none of the 20 scanned (α, β) combinations moves Ω_m toward the data-preferred 0.347, and the stated local analytic approximation explains the rigidity; the full (α, β) space has not been exhausted. Within this restricted scan, the remaining dials are (θ, m_FP, Ω_Λ), and the tested variations incur compressed-CMB penalties: raising θ to reach Ω_m = 0.315 costs Δχ²_CMB ≈ 31.65 relative to the free P1 ΛCDM reference. Varying the anchor Ω_Λ fares no better: along the degeneracy ray (closure Ω_m·cos²θ = 1 − Ω_Λ − Ω_r), at Ω_Λ = 0.657 the model reaches the SNe best-fit region (Ω_m = 0.347, Δχ²_SNe = 0.05), but the CMB Δχ² rises by ≈ 30.0 and the BAO Δχ² by ≈ 24.4 (Δχ²_BAO(Ω_Λ) has a sharp minimum at the baseline anchor: Δχ²_BAO = 24.61 at Ω_Λ = 0.657 and 21.99 at 0.740, vs 0.20 at the baseline 0.699, block-diagonal scheme). The combined cost of moving along the ray is ≈ 47.9: the anchor is a sharp minimum along the tested closure ray under the stated per-probe references, not a hiding spot within that scan. Second, the SNe tension persists at Δχ² ≈ 13–16 across the entire (α, β) grid. The Ω_m tension is therefore structural in the restricted sense of the tested grid, baseline anchor and closure ray. A global claim over the full bimetric parameter space is not established; a full exploration of (θ, m_FP, Ω_Λ) remains necessary.
The self-consistent Ω_m = 0.304 differs from both the Planck value (0.315) and the DESI BAO best-fit (0.2972). This is not a discrepancy but a prediction: the quartic constraint determines Ω_m from the bimetric parameters, not from data. The H₀ = 68.84 km/s/Mpc follows from the BAO distance fit with f = 29.608 and Planck r_d = 147.09 Mpc (Planck 2018 TT,TE,EE+lowE best-fit [Planck2018], consistent with the DESI DR2 fiducial value r_d = 147.05 Mpc adjusted for ω_b = 0.02237), consistent with the joint analysis of Högås & Mörtsell (2025), who obtained H₀ = 69.0 ± 0.4 km/s/Mpc. The bimetric model does not claim to resolve the Hubble tension in this analysis: Lüben et al. (2020) [Luben2020] showed that resolution requires combining CMB, SNe Ia, and BAO data, where the phantom phase in the z ~1–10 interval shifts the inferred H₀ toward the local value. We note that Högås & Mörtsell (2026) [Hogas2026] have shown that physically motivated priors in the local distance ladder reduce the Hubble tension from 5σ to 2σ within ΛCDM.
To quantify whether the Ω_m tension can be alleviated within the bimetric parameter space, we scan Ω_m as a function of the mixing angle θ (fixing α = 0, β = 3, m_FP/H₀ = 25). Increasing θ from 6.1° to 20° raises the predicted Ω_m from 0.304 to 0.341: the larger mixing angle shifts y(0) further from unity, which increases the subtracted massive-graviton term and thereby decreases the effective dark-energy density, raising Ω_m = 1 - Ω_r - Ω_DE(0). However, this alleviation comes at a steep cost: the absolute and comparative CMB χ² grows from 0.041 to 370.45 over the same range (free P1 reference) (Table A5). The Planck value Ω_m = 0.315 is reached at θ ≈ 12.2°, where Δχ²_CMB ≈ 31.65 within the stated geometric convention, and the Pantheon+ full-covariance value Ω_m = 0.334 requires θ ≈ 18.3°, where Δχ²_CMB ≈ 244.80. Within the tested conventions, reconciling Ω_m with these external estimates requires mixing angles that are strongly penalised by the compressed-CMB geometric diagnostic. This is not a global exclusion by a full bimetric CMB calculation.
### 7.6. BIC Comparison
The Bayesian Information Criterion for BAO gives ΔBIC_min = 0.188 (from the (α, β) scan minimum, matching but not derived from the baseline Δχ² = 0.188) (bimetric parameters fixed, δ k = 0) — models are indistinguishable for BAO. For the hypothetical full-model count (δ k = 4), ΔBIC_full = 0.188 + 4ln(13) = 10.4 in this descriptive BAO-only calculation. This BIC applies to BAO only and is neither evidence nor a Bayes factor; it is not used below as a joint model-selection result. In the BAO analysis, all four bimetric parameters are held fixed; the BIC with δk = 4 is a hypothetical estimate for the full model, not a fit to BAO data. For SNe, Δχ² = 15.57 (6.50 with the full covariance) disfavours the bimetric model at the baseline point, though the diagonal-covariance value is partly driven by the diagonal covariance approximation. At the CMB baseline, Δχ² = 0.041 against the free P1 reference (χ²_ref ≈ 0): it adds only a small positive penalty, with no standalone model-selection significance assigned.
### 7.7. Forecasts: Late-ISW and S₈
Two conditional research directions follow from the perturbatively stable window z < z* (= 7.880536). Late-ISW: the unscreened regime below z* sources an Integrated Sachs–Wolfe signal accumulated along the line of sight; the forecast is conditional on the bimetric Boltzmann treatment in the surviving pocket (§6.4). S₈: the growth suppression implied by the Ω_m prediction (0.304 vs 0.347) propagates to the amplitude of structure; the forecast is likewise conditional. Both forecasts are stated for the stable window only — no claim is made for z > z*.
**Background-level directions only.** The effective background motivates further study of late-ISW and S₈, but this draft does not assign quantitative shifts or observational preferences to them. Only an appropriate perturbative/Boltzmann computation can turn these directions into predictions.
### 7.8. Limitations and Future Work
This public draft has the following explicit limits:

1. **Restricted per-probe scope.** BAO blindness is statistical indistinguishability at the baseline and across the reported grids, under the stated covariance and profiling conventions. It is not a theorem for the entire bimetric theory or for future data.
2. **Profiling is not marginalisation.** The scale f and amplitude M are profiled by minimisation. The quantitative effect of this choice, the finite interval 20 ≤ f ≤ 80, and boundary hits on the scientific conclusions remain under audit. No evidence or Bayes factor is calculated.
3. **SNe covariance.** The full Pantheon+ STAT+SYS result Δχ² = 6.50 is primary. The diagonal value 15.57 and the σ_int scan are sensitivity diagnostics; σ_int is not added on top of the full covariance.
4. **GR-calibrated compressed CMB.** The compressed likelihood supplies a conditional geometric background diagnostic with fixed z* = 1090, not a full bimetric Boltzmann, ISW, lensing, or perturbation test. The physical error of this transfer is unknown. The earlier numerical claim that perturbative corrections must shift Δχ²_CMB by ≲ 10⁻³ is withdrawn: the available evidence does not establish that bound.
5. **Legacy neutrino conventions.** P1/P2 and their piecewise neutrino treatments belong to the reported legacy scans. The developing P3 exact-Fermi–Dirac background is a different target and is not used to claim updated numerical results here.
6. **No accepted joint P3 result.** Exploratory historical posteriors are not accepted P3 posteriors. There is no accepted production MCMC, full G3, joint parameter constraint, matched model-selection procedure, evidence, or Bayes factor in this draft.
7. **Computational domain and control.** The effects of root/closure rules, numerical unknowns, prior coordinates, and matched ΛCDM control remain under methodological audit. Failure of an optimiser or a numerical closure is not by itself a physical falsification.
8. **Baseline scalar instability.** The gradient-instability and growth statements apply to the adopted finite branch, baseline parameters, and linear/WKB diagnostic. They do not determine nonlinear saturation or viability across the full parameter space.
9. **Version-specific reproducibility.** Executable code v10.6 reproduces the legacy v21.5/v21.6 pipeline, not every updated v28.0 value. The documentary archive v13.3 and a mapping of evidence are supplied separately; no claim of one-command reproduction of the entire manuscript is made.
10. **Conditional forecasts.** Future-survey, late-ISW, S₈, and H₀ statements are possibilities or background extrapolations, not certified predictions. No CMBComp adoption or bimetric Boltzmann calculation is included in this release.
## 8. Conclusions
This public working draft presents a restricted per-probe decomposition, not a final joint comparison of bimetric gravity and ΛCDM.

1. **BAO.** At the baseline point and across the reported grids, DESI DR2 BAO differences are small: Δχ² = 0.188 with the full covariance at baseline. “Blindness” is limited to this statistical and computational domain.
2. **Mechanism.** The late-time bimetric correction is small, and profiling of f absorbs part of the residual scale dependence. The pure-correction and Ω_m-mismatch diagnostics separate these contributions; profiling is not Bayesian marginalisation.
3. **SNe.** Pantheon+ gives Δχ² = 6.50 with full STAT+SYS covariance and 15.57 with diagonal covariance. These fixed-benchmark, per-probe differences do not constitute a joint model-selection significance.
4. **Ω_m rigidity.** The near-rigidity of Ω_m is established for the reported 20-point (α, β) grid and stated local analytic approximation. A global impossibility claim over the full interaction space is not made.
5. **Compressed CMB.** Within the GR-calibrated geometric convention, the θ penalty rises from 0.04097 at baseline to 2369.29001 at 30° in P1. This is not a full bimetric CMB exclusion, and the physical error of the approximation is unquantified.
6. **Effective background.** The benchmark exhibits an effective phantom w_DE in the stated redshift regime. Quantitative consequences for H₀, S₈, or ISW require a dynamical calculation and are not claimed here.
7. **Scalar sector.** At the baseline and on the adopted finite branch, the linear/WKB diagnostic gives z* = 7.880536 and very large growth over the unstable window. This does not decide viability at other parameter points or after nonlinear evolution.
8. **Joint analysis.** The exact-Fermi–Dirac P3 program remains under audit. No accepted production MCMC, full G3, matched model-selection procedure, evidence, or Bayes factor is presented. Earlier public code v10.6 remains available as a legacy reproducible pipeline; version-specific coverage is documented in the release package.
## Acknowledgments
The computational work in this study was performed with the assistance of AI tools. The early main numerical computations were carried out with Alice AI (Yandex). The current phase — administration of the project workflow, code development, and numerical computations, including the quartic solver, χ² minimisation, B-parameter evaluation, and the full-covariance supernova pipeline — was carried out using Python code developed in collaboration with Kimi (Moonshot AI), in the role of the project's AI administrator ("Алиса"). The manuscript underwent multiple rounds of review by AI-assisted auditors (GigaChat Standard, GigaChat Reasoning, and GigaChat Research, Sber), whose feedback identified formatting issues, logical inconsistencies, and bibliographic errors. Independent computational cross-checks of key results were performed by additional AI systems used as verifiers: DeepSeek (DeepSeek AI) confirmed all seven verified quantities (B-parameters, quartic root, Ω_m self-consistency, H₀, ΔBIC, and Δχ² decompositions) to full precision, and — decisively for the CMB results — the compressed-likelihood θ-scan was independently reproduced by three separate, self-contained codebases written by DeepSeek, Qwen (Alibaba), and GLM (Zhipu AI), each executed verbatim by the project administrator; all three reproduced the canonical values within 0.1% (with one of these cross-checks localising a real sign error in the sound-speed formula, fixed in archive v10.6). The perturbative-stability result of §6 (the mapping α_A = tan θ, m²β_n^(A) = B_n, and the X(z) table) was verified in a blind series: independent code written by DeepSeek (DeepSeek AI), GLM (Zhipu AI), and Kimi (Moonshot AI) — none of which had access to our results or to each other's work — reproduced the sign pattern, the crossing redshift z* ≈ 7.9, and X(1090) ≈ −1.16 to within their stated numerical conventions; a fourth participant (Qwen, Alibaba) initially failed verification (its reported table was not produced by its submitted code) and reproduced the result only after a partial hint, and is therefore not counted as independent. One of the blind responses (DeepSeek) identified a genuine error in our initial mapping (the consistency condition m²β_n = B_n), which was corrected and re-verified before inclusion. The full-covariance supernova analysis was independently reproduced from the public Pantheon+ data before inclusion (§4.2). AI-generated numerical content used in this paper was subjected to the project's documented checks; in particular, numerical claims supplied by external AI systems were accepted only after verbatim execution of their code, never on the basis of reported values alone. The author takes full responsibility for all results and conclusions presented in this work.
## Data Availability
The DESI DR2 BAO data are publicly available at https://data.desi.lbl.gov/doc/releases/dr2/ [DESI2025a, DESI2025b]. The Pantheon+ data are publicly available at https://github.com/PantheonPlusSH0ES/DataRelease [Brout2022, Scolnic2022]. The y(z), B-parameter, w_DE(z), and scan tables are provided in this paper.

The public repository is https://github.com/chernicovalex-prog/bao-blindness-bimetric. Executable archive v10.6, accompanying manuscripts v21.5/v21.6, remains available with DOI 10.5281/zenodo.23292160. Manuscript v28.0 is a public working draft with version-specific reproducibility documented in `REPRODUCIBILITY.md`; documentary code archive v13.3 is included without claiming that the legacy executable reproduces all updated values. External cross-checks are accepted only after verbatim code execution, with logs and hashes retained. Pantheon+ data are not redistributed in the repository.
## References
[Hassan2012] S. F. Hassan and R. A. Rosen, "Bimetric Gravity from Ghost-Free Massive Gravity," JHEP 02, 126 (2012), arXiv:1109.3515.
[vonStrauss2012] M. von Strauss, A. Schmidt-May, J. Enander, E. Mörtsell, and S. F. Hassan, "Cosmological Solutions in Bimetric Gravity and their Observational Tests," JCAP 03, 042 (2012), arXiv:1111.1655.
[Comelli2012] D. Comelli, M. Crisostomi, F. Nesti, and L. Pilo, "FRW Cosmology in Ghost Free Massive Gravity," JHEP 03, 067 (2012), Erratum: JHEP 06, 020 (2012), arXiv:1111.1983.
[Kunz2009] M. Kunz, "The dark degeneracy: On the number and nature of dark components," Phys. Rev. D 80, 123001 (2009), arXiv:astro-ph/0702615.
[Hogas2021] F. Högås and E. Mörtsell, "Analytical constraints on bimetric gravity," JCAP 05, 001 (2021), arXiv:2101.08794.
[Hogas2021b] F. Högås and E. Mörtsell, "Observational Constraints on Bimetric Gravity," JCAP 05, 002 (2021), arXiv:2101.08795.
[Hogas2025] F. Högås and E. Mörtsell, "Bimetric gravity improves the fit to DESI BAO and eases the Hubble tension," Phys. Rev. D 112, 103515 (2025), arXiv:2507.03743.
[Luben2020] M. Lüben, A. Schmidt-May, and J. Smirnov, "Vainshtein Screening in Bimetric Cosmology," Phys. Rev. D 102, 123529 (2020), arXiv:1912.09449.
[Luben2020b] M. Lüben, A. Schmidt-May, and J. Weller, "Physical parameter space of bimetric theory and SNIa constraints," JCAP 09, 024 (2020), arXiv:2003.03382.
[Konnig2014] F. Könnig, Y. Akrami, L. Amendola, M. Motta, and A. R. Solomon, "Stable and unstable cosmological models in bimetric massive gravity," Phys. Rev. D 90, 124014 (2014), arXiv:1407.4331.
[Konnig2015] F. Könnig, "Higuchi Ghosts and Gradient Instabilities in Bimetric Gravity," Phys. Rev. D 91, 104019 (2015), arXiv:1503.07436.
[Bassi2023] A. Bassi, M. Adil, M. Rajvanshi, and A. A. Sen, "Cosmological Evolution in Bimetric Gravity," arXiv:2301.11000 (2023).
[Hogas2019stab] M. Högås, F. Torsello, and E. Mörtsell, "On the stability of bimetric structure formation," JCAP 03, 021 (2020), arXiv:1910.01651.
[DESI2025a] DESI Collaboration, "DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints," Phys. Rev. D 112, 083515 (2025), arXiv:2503.14738.
[DESI2025b] DESI Collaboration, "DESI DR2 Results I: Baryon Acoustic Oscillations from the Lyman Alpha Forest," Phys. Rev. D 112, 083514 (2025), arXiv:2503.14739.
[Eisenstein2005] D. J. Eisenstein et al., "Detection of the Baryon Acoustic Peak in the Large-Scale Correlation Function of SDSS Luminous Red Galaxies," Astrophys. J. 633, 560 (2005), arXiv:astro-ph/0501171.
[Brout2022] D. Brout et al., "The Pantheon+ Analysis: Cosmological Constraints," Astrophys. J. 938, 110 (2022), arXiv:2202.04077.
[Scolnic2022] D. Scolnic et al., "The Pantheon+ Analysis: The Full Data Set and Light-curve Release," Astrophys. J. 938, 113 (2022), arXiv:2112.03863.
[Bansal2025] P. Bansal and D. Huterer, "Expansion-history preferences of DESI DR2 and external data," Phys. Rev. D 112, 023528 (2025), arXiv:2502.07185.
[Akrami2015] Y. Akrami, S. F. Hassan, F. Könnig, A. Schmidt-May, and A. R. Solomon, "Bimetric gravity is cosmologically viable," Phys. Lett. B 748, 37 (2015), arXiv:1503.07521.
[HuSugiyama1995] W. Hu and N. Sugiyama, "Anisotropies in the Cosmic Microwave Background: An Analytic Approach," Astrophys. J. 444, 489 (1995), arXiv:astro-ph/9407093.
[Hogas2026] F. Högås and E. Mörtsell, "Physically-motivated priors in the local distance ladder significantly reduce the Hubble tension," MNRAS 548, stag724 (2026), arXiv:2601.22215.
[Planck2015XIV] Planck Collaboration, "Planck 2015 results. XIV. Dark energy and modified gravity," Astron. Astrophys. 594, A14 (2016), arXiv:1502.01590.
[Planck2018] Planck Collaboration, "Planck 2018 results. VI. Cosmological Parameters," Astron. Astrophys. 641, A6 (2020), arXiv:1807.06209.
## Appendix A: Supplementary Tables
### Table A1: B-Parameters
Baseline parameters: θ = 6.1°, m_FP/H₀ = 25, α = 0, β = 3, Ω_Λ = 0.699, Ω_r = 9.08 × 10⁻⁵.
```
Parameter   Value       Formula
B0          -40.248     3ΩΛ − sin²θ·m²_FP·(3 + 3α + β)
B1          28.230      sin²θ·m²_FP·(1 + 2α + β)
B2          -21.173     −sin²θ·m²_FP·(α + β)
B3          21.173      sin²θ·m²_FP·β
B4          -28.206     3tan²θ·ΩΛ + sin²θ·m²_FP·(−1 + α − β)
B1²         796.94      B1² = 796.94 (full-precision); 796.93 from displayed B1
```
### Note on the Ω_m conventions (0.30435 / 0.3029 / 0.30312)
The predicted Ω_m appears in three numerically close forms, all reproduced by direct re-runs of the pipelines (reference closure, archive v13.2); the differences are convention, not bug. The closure relation is Ω_m = 1 − Ω_r − ω_ν/h² − Ω_DE(0), with the neutrino term ω_ν/h² (Σm_ν = 0.06 eV) the only moving part: (i) Ω_m = 0.30435 in the BAO/SNe pipeline, where no explicit ω_ν term enters closure (neutrinos enter only through Ω_r in the dynamics); (ii) Ω_m = 0.3029 when the ω_ν/h² term is included and evaluated at the reference h = 0.678 (Table A1b below); (iii) Ω_m = 0.30312 in the CMB θ-scan table (archive v10.6, table 9), where the same term is evaluated at the fitted h = 0.72425 of the scan iteration. All three recover E²(0) = 1 to ~10⁻¹⁵; the spread (≲ 1.5 × 10⁻³) is the neutrino closure contribution at two different h.

### Table A1b: Rigidity of the Ω_m Prediction over the (α, β) Grid

Ω_m predicted by the quartic self-consistency condition at z = 0, θ = 6.1°, for the 20-combination (α, β) scan. The prediction is independent of (α, β) to four decimal places: the α, β-dependent terms in Ω_DE(0) are of order (1−y(0))² ≈ 2 × 10⁻⁷. Consequently, no choice of potential parameters within this class can reconcile the predicted Ω_m = 0.304 with the SNe-preferred 0.347 (full covariance) — the Ω_m tension is structural, not parametric.
```
alpha \ beta     1        2        3        4
-1.0           0.3029   0.3029   0.3029   0.3029
-0.5           0.3029   0.3029   0.3029   0.3029
 0.0           0.3029   0.3029   0.3029   0.3029
+0.5           0.3029   0.3029   0.3029   0.3029
+1.0           0.3029   0.3029   0.3029   0.3029
```

### Table A2: y(z) and E²(z) — Bimetric vs ΛCDM
Fractional correction δE² = (E²_bim − E²_LCDM)/E²_LCDM. Maximum |δE²| ≈ 1.11% at z ≈ 5. The correction is negative (bimetric expansion is slower) at all z > 0.
```
z       y(z)      E2_bimetric    E2_LCDM     dE2 (%)
0.00    0.999512  1.000000       1.000000    -0.0000
0.10    0.999350  1.099643       1.100784    -0.1036
0.50    0.998352  1.715027       1.723211    -0.4750
1.00    0.996093  3.107718       3.131845    -0.7704
2.00    0.986814  8.830855       8.920487    -1.0048
2.33    0.981968  11.821478      11.945333   -1.0368
3.00    0.968770  19.980325      20.197501   -1.0753
5.00    0.896071  65.815888      66.553851   -1.1088
10.00   0.556527  403.379176     407.121089  -0.9191
12.00   0.442652  666.804264     671.956212  -0.7667
1090    0.000002  523876594.328  523876608.439 -0.0000
```
### Table A3: Effective Dark Energy Equation of State w_DE(z)
Phantom regime (w_DE < −1) at all z where Ω_DE is significant (z ≲ 5). Ω_DE → 0 at z ≈ 4.9, making w_DE ill-defined beyond that point.
```
z       w_DE(z)
0.0     -1.0050
0.5     -1.0169
1.0     -1.0411
2.0     -1.1536
3.0     -1.4610
4.0     -2.5988
4.5     -5.4906 † (ΩDE -> 0; ill-defined)
5.0     16.2874 † (ΩDE -> 0; ill-defined)
```
† Beyond the zero crossing of Ω_DE(z) at z ≈ 4.9 the ratio w_DE = p_DE/ρ_DE is a ratio of two vanishing quantities and carries no physical meaning; the negative-Ω_DE asymptote is discussed in §7.2.
### Table A4: D_M(z_*) Accumulation by Redshift
Fraction of the total comoving distance D_M(z_* = 1090) contributed by z ≤ z_check. The unscreened region (z ≲ 12, where H < m_FP) contributes 72.19% of D_M(z_*), despite spanning only ~1% of the redshift range to z*.
```
z <= z_check    % of D_M(z_*)    Regime
1               24.37%          Unscreened
5               57.20%          Unscreened
10              69.43%          Unscreened
12              72.19%          Screening boundary
50              87.85%          Screened
500             98.51%          Screened
1090            100.00%         z_* (last scattering)
```
### Table A5: CMB θ-Scan — Absolute χ² and Δχ² vs Mixing Angle
CMB-only values in the unified independently reproduced vintage. The free (ω_m, h, ω_b) ΛCDM reference has χ²_ref ≈ 0 and DOF = 0 in each convention; therefore absolute and comparative values coincide to the displayed precision. Bimetric absolute p-values at baseline use DOF = 1; no reference goodness-of-fit p-value is defined.
```
θ (deg)   P1 χ²_bim   P1 Δχ²      P2 χ²_bim   P2 Δχ²
6.1       0.04097     0.04097     0.07712     0.07712
10.0      9.35910     9.35910     8.90408     8.90408
12.2      31.64872    31.64872    30.80732    30.80732
15.0      93.70686    93.70686    92.25547    92.25547
18.3      244.79718   244.79718   242.44979   242.44979
20.0      370.44548   370.44548   367.55851   367.55851
25.0      1029.62953  1029.62953  1024.83014  1024.83014
30.0      2369.29001  2369.29001  2362.06129  2362.06129
```
At θ = 30° with m_FP/H₀ = 25, the static Higuchi bound is satisfied (625 > 1.4). The θ-scan is a sensitivity map, not a set of independent look-elsewhere-corrected hypothesis tests.
### Table A6: SNe σ_int Scan — Δχ² vs Intrinsic Dispersion (diagonal covariance)
Diagonal covariance with σ_int added in quadrature: σ²_total = σ²_diag + σ²_int. Ω_m^best refers to the ΛCDM fit.
```
σ_int (mag)   Ω_m^best (ΛCDM)   Δχ² (standard)   Δχ² (pure)   χ²/dof
0.00          0.3742             15.57            1.49         0.455
0.05          0.3743             14.68            1.41         0.428
0.08          0.3743             13.49            1.29         0.393
0.10          0.3743             12.55            1.20         0.365
0.12          0.3742             11.58            1.11         0.338
0.15          0.3740             10.15            0.98         0.298
0.20          0.3735             8.06             0.78         0.239
```

## Appendix B: Pseudocode for E²(0) Residual Computation
The E²(0) residual (~10⁻¹⁵) serves as a numerical consistency check. The following pseudocode describes the computation:
```
Input: θ, m_FP, α, β, Ω_Λ, Ω_r
1. Compute B-parameters:
   B₀ = 3Ω_Λ − sin²θ·m²_FP·(3 + 3α + β)
   B₁ = sin²θ·m²_FP·(1 + 2α + β)
   B₂ = −sin²θ·m²_FP·(α + β)
   B₃ = sin²θ·m²_FP·β
   B₄ = 3tan²θ·Ω_Λ + sin²θ·m²_FP·(−1 + α − β)
2. Find self-consistent Ω_m:
   For trial Ω_m in [0.01, 0.99]:
     a. At z = 0: Ω_tot = Ω_m + Ω_r
     b. Build quartic coefficients a₀..a₄ from B-params
     c. Solve a₀ + a₁y + a₂y² + a₃y³ + a₄y⁴ = 0 using Brent's method on [10⁻⁶, 1.0]
     d. Compute E²(0) = Ω_m + Ω_r + Ω_DE(y(0)),
        Ω_DE = Ω_Λ − sin²θ·m²_FP·(1−y)·[1 + α(1−y) + β/3·(1−y)²]
   Find Ω_m such that E²(0) = 1 (Brent's method)
3. Compute residual: E²(0)_residual = |E²(0) − 1.0|; expected ~1e-15
4. Round-trip test: (θ, m_FP, α, β, Ω_Λ) -> B -> y -> recover params; ~1e-13
```

## Appendix B.2: Analytic Late-Time Expansion and the Rigidity of Ω_m

At z = 0 the de Sitter attractor pins y to unity: writing y = 1 − ε, the quartic
F(ε) ≡ f(1−ε; Ω_m(ε)) with the closure Ω_m(ε) = 1 − Ω_r − ω_ν/h² − Ω_Λ +
sin²θ·m²_FP·ε·(1 + αε + βε²/3) expanded to first order gives

  F(0) = 1 − Ω_Λ − Ω_r − ω_ν/h²          (≈ 0.2996 at baseline)
  F'(0) = −[f'(1) − sin²θ·m²_FP]          (≈ −617 at baseline, dominated by the
                                           spin-2 mass scale m²_FP = 625)
  ε = F(0)/|F'(0)| + O(ε²)                (= 4.857 × 10⁻⁴ vs the exact
                                           4.882 × 10⁻⁴ — the linear expansion
                                           reproduces ε to 0.5%)

so ε ∝ (1 − Ω_Λ − Ω_r − ω_ν/h²)/m²_FP is set by the distance from the de Sitter
point, not by the potential shape. The (α, β)-dependence of Ω_DE(0) is
explicitly O(α·ε²) and O(β·ε³), and the implicit dependence through ε adds only
|∂Ω_m/∂α| ≈ 4 × 10⁻⁹ and |∂Ω_m/∂β| ≈ 10⁻¹² (finite differences of the exact
quartic solution at baseline). Hence the present-day matter density is a
function of (θ, m_FP, Ω_Λ) alone within this model class. Note that the rigidity
is an amplitude statement at z = 0: at z ≳ 0.5 the deviation 1−y(z) grows and
(α, β) do shape E²(z) — they are constrained through the form of the expansion
history (the pure corrections of §§3–5), not through the present-day energy budget.

## Appendix C: Reproducibility and External Verification
Numerical support combines canonical reruns, internal consistency checks, and item-specific external blind replication. External code is accepted only after verbatim execution (a report is not a run), with files, logs, and hashes retained. The list below states the actual coverage; it does not claim that every sentence or every updated number was independently reproduced by an external participant.

- **CMB free reference and split test (blind series СЗ-18).** Qwen, Alice AI and Grok independently implemented the free (ω_m, h, ω_b) ΛCDM reference, P1/P2 θ-scans and the high/low-z residual split without access to one another's work. Their code was executed verbatim; the results agree with the unified values in §5 and Appendix A5. The reference has zero residual DOF, so no reference goodness-of-fit p-value is defined.
- **Exact SNe preference map.** Both 254 × 280 grids (71,120 cells per cut) were reproduced independently from the Pantheon+ data and covariance. The largest absolute Δχ² discrepancies were 5.97 × 10⁻⁸ (z > 0.01) and 4.47 × 10⁻⁸ (z > 0.023), within printed CSV precision. Appendix D uses these verified grids, not the superseded surrogate island map.
- **Analytic closure.** An independent implementation verified the quartic derivatives, exact closure, e₁/e₂ approximations and non-circular prediction of Ω_m from (h, θ). The product C·Ω_m is an approximate finite-parameter regularity, not an invariant; the CMB observables break the apparent degeneracy.
- **Background pipeline.** DeepSeek (DeepSeek AI) confirmed seven verified quantities (B-parameters, quartic root, Ω_m self-consistency, H₀, ΔBIC, Δχ² decompositions) to full precision. The compressed-likelihood CMB θ-scan was independently reproduced by three separate codebases (DeepSeek, Qwen, GLM) within 0.1%; one of these cross-checks localised a real sign error in the sound-speed formula, fixed in archive v10.6. The CPL fits of w_DE(z) were reproduced by three external codes (Alice AI, DeepSeek, Qwen) agreeing to 10⁻³.
- **Scalar sector (blind series №8 and №10).** Participants received specifications containing no target numbers and had no access to our results or to each other's work. In series №8, GLM (Zhipu AI) reproduced the quartic root at z = 0, the full X(z) sign pattern for z ≤ 100, and the convention-mapping identity m_FP² = (1/α_A² + 1)(B₁ + 2B₂ + B₃) to 10⁻¹³, and localised the crossing at z* = 7.88053597 (Brent's method) — an 8-digit agreement with the reproducible canonical-code value 7.880536. Series №10 included failed background-level submissions and is not used as an independent confirmation in this draft; this does not alter the separately documented series №8 result.
- **On the two z* values.** Earlier versions quoted z* = 7.8836 (office pipeline); the reproducible code value is 7.880536, confirmed externally as above. The difference (4 × 10⁻³) is a grid-resolution effect in the office scan; the canonical value was unified to 7.880536 in the v26 transaction. No conclusion in this paper depends on the difference.
- **Γ-block.** The growth-integral convention (1/E, per-mode k/H₀ normalisation following Akrami et al. 2015) and the value I₁ = 0.945210110 (converged quadrature; it supersedes the grid value 0.945206871188 of rev 97–98 — archive v11, finding #106) were fixed in the project journal; the Γ table of §6.2 is pure arithmetic on I₁. An earlier 1/E² convention was withdrawn as dimensionally inconsistent and is not used.


## Appendix D: Exact SNe Reference-Sensitivity Maps

Figure D1 shows the exact separable quantity Δχ²(θ, Ω_a) = χ²_bim(θ) − χ²_ΛCDM(Ω_a), evaluated on 254 θ values × 280 reference anchors = 71,120 cells for each cut. Both panels share one colour scale and display the zero contours. These are SNe-only reference-sensitivity maps, not a joint MAP, not posterior σ-contours and not a combined BAO+CMB+SNe landscape. Negative values at deliberately non-optimal anchors mean only that the fixed bimetric point fits better than that chosen ΛCDM reference; they are not evidence for a bimetric preference.

![Figure D1](figures/FigD1_SNe_reference_sensitivity.png)

*Figure D1: Exact SNe reference-sensitivity maps for z > 0.01 and z > 0.023.*

```
Quantity                         z > 0.01             z > 0.023
N                                1588                  1367
χ²_bim(6.1°)                     1439.556548           1238.426532
ΛCDM best Ω_m                    0.346994              0.345728
ΛCDM best χ²                     1433.060252           1232.634833
baseline Δ against best          6.496297              5.791698
baseline zero anchors            0.3009099, 0.3960969  0.3009092, 0.3934589
Δ at Ω_a = 0.310 (θ range)       2.36357--2.36744      2.15926--2.16302
θ-span of χ²_bim                 0.00387512            0.00376011
```

The sign is controlled mainly by the arbitrary reference anchor, while the θ-dependence of the bimetric SNe block is very small. The zero contours are only nearly vertical, not mathematically vertical. At Ω_a = 0.310 the exact values +2.36 (z > 0.01) and +2.16 (z > 0.023) confirm the SNe component of the predicted H&M-bridge at that chosen anchor. They do not verify the full H&M joint decomposition: the CMB and BAO blocks are not recomputed here, the per-block H&M values are not published, and the earlier CMB −2.9 estimate cannot be carried into the new free-reference convention without a separate calculation. No combined map is constructed.

As an empirical diagnostic on the SNe closure track, I(θ) = cos²θ·Ω_m(θ) ranges from 0.3009179301 to 0.3011875044 (+0.089584%), while Ω_m itself grows by 36.121064%. I(θ) is not a literal term of the quartic and not an exact invariant; its near-coincidence with the lower zero anchor at small θ is not an identity. The two maps and the CMB-track regularity of §2.3 use different h(θ) tracks and conventions and are not numerically identified.

*--- End of Appendix ---*
