# BAO Blindness to Bimetric Gravity: A Per-Probe Δχ² Decomposition
Aleksandr Chernikov
*Independent Researcher*
October 2026 (v21.6)
## Abstract
We build on the dataset and model setup of Högås & Mörtsell (2025) and extend it by decomposing the constraining power of each cosmological probe individually. Using DESI DR2 BAO data (13 points), Pantheon+ SNe Ia data (1588 light curves after z > 0.01 cut), and the compressed Planck 2018 + ACT DR6 likelihood (Bansal & Huterer 2025), we test bimetric gravity against ΛCDM at the background level. BAO is effectively blind to the bimetric parameters at baseline (Δχ² = 0.178, p ≈ 0.66; here 'blindness' denotes statistical indistinguishability, Δχ² ≪ 3.84, caused by the de Sitter attractor suppressing the correction to E²(z), with f-marginalisation absorbing only the residual — see §5.1), though it becomes sensitive to θ at large angles (Δχ² → 2 at θ ≈ 31.5°; see §4.1), because the dynamical evolution drives y → 1 at late times as a de Sitter attractor, rendering the bimetric correction to E²(z) intrinsically small (≲ 1.1% throughout the BAO redshift range), and the marginalisation over f = c/(H₀ r_d) absorbs residual signal. In contrast, SNe Ia are not blind: the bimetric model predicts Ω_m = 0.304 from the quartic constraint, while the ΛCDM best-fit to Pantheon+ (diagonal covariance) gives Ω_m = 0.374, yielding Δχ²_SNe = 15.58 (diagonal covariance; an upper bound, see §6). With the full Pantheon+ STAT+SYS covariance, the ΛCDM best-fit shifts to Ω_m = 0.347 — closer to the BAO/CMB-preferred range — and Δχ²_SNe decreases to 6.50 (pure bimetric correction 0.96); the preference for ΛCDM is moderate, 2.5σ–3.9σ depending on the error model. The discrepancy arises because the negative bimetric E²(z) correction (≲ 1%) shifts the model in the wrong direction relative to the SNe-preferred Ω_m; the sensitivity is dominated by the Ω_m mismatch rather than by the E²(z) correction itself (pure bimetric Δχ² = 1.49 diagonal, 0.96 full covariance). We further show that the predicted Ω_m is rigid against the potential parameters: it is independent of (α, β) to |∂Ω_m/∂α| ≲ 5 × 10⁻⁹ across the scanned grid, so the Ω_m tension is structural at the fixed joint-fit anchor, not parametric. Varying the anchor Ω_Λ itself only trades the SNe tension for CMB and BAO tension: at Ω_Λ = 0.657 the SNe Δχ² drops from 6.50 to 0.05, but the CMB Δχ² rises by ≈ 28.6 and the BAO Δχ² by ≈ 23.4 — the degeneracy ray costs ≈ 45.6 overall and offers no hiding spot (§5.5). CMB — within the GR-calibrated compressed distance prior, used here strictly as a geometric background diagnostic (§3.3; a formal CMB constraint requires a bimetric-modified Boltzmann calculation) — is also sensitive to the mixing angle θ (Δχ² grows from 0.04 to 2369 as θ increases from 6.1° to 30°, v10.6), because the comoving distance D_M(z_*) integrates the bimetric correction from z = 0 to z_* = 1090, accumulating the unscreened late-universe contribution (z ≲ 12, where H < m_FP, contributing 72% of D_M(z_*) despite spanning only 1% of the redshift range to z_*). The CMB-only sensitivity to θ, isolated from BAO and SNe contributions, is a new result of this work. Depending on the choice of error model, the combined Δχ² across the three probes is 6.7 (all-full-covariance scheme) to 16.2 (all-diagonal scheme). We quantify the blindness mechanism for each probe, provide B-parameter values (B_0 = -40.248, B_1 = 28.230), verify B-parameter round-trip consistency to 10⁻¹³ (full-precision internal computation; the displayed rounded values recover to ~10⁻⁴) and E²(0) to 10^-15, and discuss the applicability of compressed CMB likelihoods to modified gravity. We note that the small bimetric correction (≲ 1.1%) at BAO redshifts, combined with the marginalisation over the distance ratio f = c/(H₀ r_d), makes θ-constraint through BAO challenging at current precision, though not impossible in principle for future sub-percent surveys.
**Keywords: bimetric gravity, DESI BAO, dark energy, modified gravity, CMB, supernovae**
## 1. Introduction
Bimetric gravity (Hassan & Rosen 2012 [Hassan2012]) extends general relativity by introducing a dynamical second metric f_μν coupled to the physical metric g_μν through a Fierz--Pauli mass term. The theory is ghost-free at the fully nonlinear level and admits viable cosmological solutions on an FLRW background (von Strauss et al. 2012 [vonStrauss2012]; Comelli et al. 2012 [Comelli2012]).
Högås & Mörtsell (2021) [Hogas2021] derived the B-parameter formalism (observational constraints in the companion paper [Hogas2021b]) for bimetric cosmology, and Högås & Mörtsell (2025) [Hogas2025] provided the full Friedmann equation, Ω_DE(z), and quartic solver for the scale factor ratio y(z) = a_f/a_g (following the Appendix A definition of Högås & Mörtsell (2025), Eq. A10; note that the main text of Högås & Mörtsell 2025 contains a typographical error where y is defined as a_g/a_f, which is inconsistent with their own Eq. A10 and with the physical requirement y → 0 at the Big Bang). Their analysis combined BAO, CMB, and SNe Ia data through MCMC to obtain joint constraints on the bimetric parameters. However, they did not decompose the constraining power of individual probes---i.e., they did not identify which probe drives the constraints and which are blind.
This work builds on the dataset and model setup of Högås & Mörtsell (2025) and extends it by decomposing the constraining power of each probe individually. Table 1 summarises the key differences. Our original contributions are: (1) the BAO χ² pipeline with marginalisation over r_d, (2) numerical B-parameter values, (3) the probe-by-probe Δχ² decomposition revealing that BAO is blind while SNe and CMB are sensitive, (4) the χ² decomposition for BAO, (5) the round-trip consistency test, (6) the CMB-only θ sensitivity as an isolated result, and (7) the full-covariance SNe assessment (Δχ² = 6.50), which closes the leading quantitative limitation of the diagonal analysis.
Lüben, Schmidt-May & Smirnov (2020) [Luben2020] demonstrated that the early universe in bimetric theory is screened by the Vainshtein mechanism on an FLRW background, with the spin-2 mass setting the cosmological Vainshtein scale. In a global analysis using CMB, Cepheid, SNe Ia, and BAO data, they showed that a minimal β_0β_1β_4 bimetric model can resolve the H₀ tension via an effective phantom dark energy component (w_DE < -1) in the redshift interval 1 ≲ z ≲ 10. Their finding that current BAO data do not distinguish bimetric models from ΛCDM is consistent with and complementary to our analysis: we quantify this indistinguishability per probe and identify its origin in the smallness of the bimetric correction at late times.
```
Table 1. Comparison with Högås & Mörtsell (2025).

Aspect                  H&M 2025                This work
Data                    BAO + CMB + SNe Ia      Each probe individually
                        (combined)              
Method                  MCMC (joint)            Grid scan + Brent (per probe)
Omega_m                 Derived from quartic    Predicted by quartic
                        (varies in MCMC)        
B-parameters            Not reported            Numerical (§2.4)
                                                numerically             
Delta chi2              Not performed           Per probe (Table 2)
decomposition                                   
BAO blindness           Not quantified          Delta chi2 = 0.178 (full, incl.
                                                Omega_m penalty; pure E2 corr.
                                                0.023), flat across
                                                20 (alpha,beta) combos and 4
                                                decades in m_FP/H0 (0.023 pure
                                                bimetric, see §4.1)
SNe sensitivity         Not quantified          Delta chi2 = 15.58 (diag), 6.50
                                                (full cov); 1.49/0.96 (pure);
                                                11.58 with sigma_int = 0.12
CMB sensitivity         Not isolated            Delta chi2 = 0.04--2369
to theta                                        
theta sensitivity       Not provided            grows to ~2 at theta ~ 31.5 deg
(BAO only)                                      (block-diagonal covariance)
Round-trip test         Not performed           10^-13 (B-param), 10^-15 (E2)
```
## 2. Theoretical Framework
### 2.1. Bimetric Friedmann Equation
The bimetric Friedmann equation for the physical metric g_μν takes the form (Högås & Mörtsell 2025, Eq. A12): E²(z) ≡ H²(z)/H₀² = Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_DE(z), where Ω_DE(z) = Ω_Λ - sin²θ · m_FP² · (1-y)[1 + α(1-y) + β/3(1-y)²], where θ is the mixing angle, m_FP is the Fierz--Pauli mass, α and β are potential parameters, and Ω_Λ is the effective cosmological constant, treated as a free model parameter (following Högås & Mörtsell 2025, who parameterize it via ω_Λ = Ω_Λ h²). At the baseline parameters (α = 0, β = 3, m_FP/H₀ = 25, θ = 6.1°), Ω_Λ ≈ 0.699. In the GR limit (y → 1, θ → 0), Ω_DE → Ω_Λ and the standard Friedmann equation is recovered. The physical interpretation is that Ω_DE receives a negative contribution proportional to sin²θ · m_FP², representing the energy density of the massive graviton. As y deviates from 1 (the GR limit), this correction grows, modifying the expansion history.
### 2.2. Quartic Equation for y(z)
The scale factor ratio y(z) satisfies a quartic equation (Högås & Mörtsell 2025, A13; corresponding to Eq. 5.5 of Högås & Mörtsell 2021 [Hogas2021]): a_0 + a_1 y + a_2 y² + a_3 y³ + a_4 y⁴ = 0, with coefficients per Eq. A13 of Högås & Mörtsell (2025), where Ω_tot(z) = Ω_m(1+z)³ + Ω_r(1+z)⁴ and m_FP² is in units of H₀². The coefficients follow the B-parameter formalism of Högås & Mörtsell (2021, 2025), solving for y = a_f/a_g directly. We solve numerically using Brent's method on the interval [1e-6, 1.0] (the lower end chosen to avoid the singular limit while containing the full physical branch), which robustly finds the physical (finite) branch: y → 1 as z → 0 (GR limit) and y → 0 as z → ∞ (early universe), consistent with the finite branch definition of Högås & Mörtsell (2021, Sec. 5.1) [Hogas2021]. The infinite branch, where y grows with z, is excluded due to the Higuchi ghost (Högås & Mörtsell 2021, Sec. 5.1). The solution is verified by computing E²(0) = 1.000000000000005, with a residual of ~5 × 10⁻¹⁵ (consistent with the self-consistency condition of Högås & Mörtsell 2021 [Hogas2021], Eq. 5.6). Note that E²(0) = 1 exactly by the self-consistency condition Ω_m = 1 - Ω_r - Ω_DE(0); the residual ~5 × 10⁻¹⁵ reflects the numerical precision of Brent's method in solving the quartic equation and serves as a numerical consistency check.
*Figure 1: Fractional bimetric correction to the expansion history, E²(z)/E²_ΛCDM*
### 2.3. Self-Consistent Ω_m
The bimetric model predicts Ω_m from the quartic constraint, rather than fitting it as a free parameter. At the baseline parameters (α = 0, β = 3, m_FP/H₀ = 25, θ = 6.1°), the quartic yields y(0) = 0.99951 and Ω_m = 0.30439. The fractional bimetric correction to E²(z) is shown in Fig. 1 and tabulated in Table A2.
### 2.4. B-Parameters
The B-parameter formalism (Högås & Mörtsell 2021 [Hogas2021]) provides an alternative parametrisation of the bimetric model, where the interaction potential parameters β_n (n = 0,...,4) are re-expressed through observables B_0 and B_1 with direct physical interpretations: B_0 characterises the effective cosmological constant and B_1 the massive graviton contribution. The parameters α and β used throughout this work are reparametrised combinations of the β_n (see Högås & Mörtsell 2021 for the exact mapping). The numerical values (Table A1) are computed from the model parameters (θ, m_FP, α, β, Ω_Λ) via the B-parameter definitions in Högås & Mörtsell (2021, Eq. 3.4) [Hogas2021], which express B_0 and B_1 as functions of the interaction potential parameters β_n and the cosmological parameters. B_0 = -40.248, B_1 = 28.230, B_2 = -21.173, B_3 = 21.173, B_4 = -28.206, B_1² = 796.94 (from full-precision value; 28.230² = 796.93 at displayed precision). All B-parameter values are recovered to a precision of 10⁻¹³ in the internal computation, while E²(0) = 1.000000000000005 is recovered to 10⁻¹⁵, confirming the consistency of the formalism. The round-trip test verifies that the full computational pipeline is self-consistent: starting from the model parameters (θ, m_FP, α, β, Ω_Λ), we solve the quartic for y(z), compute E²(z), reconstruct the quartic coefficients, and recover the original model parameters to a precision of 10⁻¹³. This ensures that no information is lost or corrupted in the forward--backward computation cycle. Note: the displayed B-parameter values are rounded to 3--6 significant figures; exact recovery from these rounded values yields deviations of order 10⁻4, not 10⁻¹³. The 10⁻¹³ precision refers to the full-precision internal computation.
### 2.5. Effective Equation of State
The effective dark energy equation of state is w_DE(z) < -1 (phantom) at all redshifts where Ω_DE is significant (i.e., z ≲ 5; at z ≈ 5, Ω_DE → 0 and w_DE becomes ill-defined), with w_DE(0) = -1.005 (Table A3). This is a genuine prediction of the bimetric model, arising from the y-dependence of Ω_DE. In the bimetric framework, this phantom regime arises naturally from the massive graviton dynamics and does not require additional phantom fields (Fig. 2).
*Figure 2: Effective dark energy equation of state w_DE*
## 3. Data and Methodology
### 3.1. DESI DR2 BAO Data
We use the DESI DR2 BAO measurements [DESI2025a, DESI2025b], building on the original BAO detection by [Eisenstein2005], consisting of 13 data points spanning redshifts z = 0.295--2.33: D_M/r_d, D_H/r_d, and D_V/r_d. The DESI DR2 covariance matrix (Cobaya release) is block-diagonal, with 2×2 blocks for each tracer pair (D_M, D_H) and a 1×1 block for D_V --- there are no cross-tracer correlations. Using exact σ values from the Cobaya matrix yields Δχ² = 0.178. Using rounded σ values from Table IV of the DESI DR2 paper yields Δχ² = 0.195 --- a 10% difference from rounding. The exact value (0.178) is adopted as the baseline. The BAO observables are computed as D_M/r_d = c/H₀ r_d ∫₀ᶻ dz'/E(z'), D_H/r_d = c/H₀ r_d · 1/E(z), D_V/r_d = [z · (D_M/r_d)² · (D_H/r_d)]^(1/3). The parameter f = c/(H₀ r_d) is marginalised over by χ² minimisation using Brent's method on the interval [20, 80], which spans the physically expected range f ≈ 25–35 (for r_d ≈ 147 Mpc and H₀ ≈ 68–74 km/s/Mpc) with generous margins.
### 3.2. Pantheon+ SNe Ia Data
We use the Pantheon+ sample (Brout et al. 2022 [Brout2022]; Scolnic et al. 2022 [Scolnic2022]), consisting of 1701 light curves from 1550 unique supernovae. After applying the standard z > 0.01 cut (removing 113 low-redshift calibrators), 1588 light curves remain in the redshift range z = 0.01--2.26. Högås & Mörtsell (2025) employ a stricter z > 0.023 cut. The choice of z > 0.01 shifts the best-fit ΛCDM Ω_m to higher values, potentially enhancing the apparent tension with the bimetric prediction Ω_m = 0.304. We assess the robustness of the SNe sensitivity to this choice via the σ_int scan (§4.2) and the full-covariance analysis (§4.2). The distance modulus is μ = 5log₁₀(D_L) + 25, where D_L = (1+z) D_M, and M is treated as a free parameter, marginalised jointly with Ω_m for the ΛCDM fit; for the bimetric model, Ω_m is predicted by the quartic constraint (§2.3) and held fixed; only M is marginalised.
Covariance treatment. We report two statistically consistent estimates. (i) Diagonal approximation (m_b_corr_err_DIAG): neglects off-diagonal systematic correlations. This approximation yields Ω_m^best = 0.374 for ΛCDM, a 2.2σ shift from the Pantheon+ full-covariance result (Ω_m = 0.334 ± 0.018 [Brout2022]), attributable to the neglected correlations between low-z and high-z supernovae. This shift directly affects Δχ²_SNe. (ii) Full covariance (STAT+SYS, 1588 × 1588 after the cut, with analytical marginalization over M): the ΛCDM best-fit moves to Ω_m = 0.347 and Δχ²_SNe decreases from 15.58 to 6.50 (§4.2). The result is robust to the low-redshift cut: with the Högås & Mörtsell cut z > 0.023 (N = 1367), the full-covariance values are Δχ²_SNe = 5.79 (pure 0.88) and Ω_m^best = 0.346. The diagonal value 15.58 should therefore be interpreted as an upper bound on the SNe tension. A σ_int scan with the diagonal covariance (Table A6) provides an intermediate, conservative estimate (11.58 at σ_int = 0.12 mag).
### 3.3. Compressed CMB Likelihood
We use the compressed Planck 2018 + ACT DR6 likelihood from Bansal & Huterer (2025) [Bansal2025], which provides a 3-parameter data vector: R = 100·√(Ω_b h² + Ω_cdm h² + Ω_ν,m h²) · D_M(z_*)/c (shift parameter, following Bansal & Huterer 2025, Eq. 4). Tensor modes in ghost-free bimetric gravity propagate at the speed of light (c_gw = c), so the model trivially satisfies the GW170817 multimessenger bound. Massive neutrinos (Σm_ν = 0.06 eV) are evolved as nonrelativistic matter at z ≲ 130 and as radiation above, following the convention of the compressed likelihood (Bansal & Huterer 2025, Appendix C); the R prefactor follows their Eq. 4 exactly, ℓ_a = π D_M(z_*)/r_s(z_*), Ω_b h² = 0.022371, with the 3×3 covariance matrix (in units of 10⁻⁸): C = [1559.83, −1325.41, −36.45; −1325.41, 714691.80, 269.77; −36.45, 269.77, 2.10] × 10⁻⁸. The sound horizon is r_s(z_*) = ∫_{z*}^∞ c_s(z)/H(z) dz, c_s(z) = c/√(3(1+R_b(z))), R_b(z) = 3Ω_b / [4Ω_γ(1+z)], where Ω_γ h² = 2.469 × 10⁻5 and z_* = 1090. The neutrino effect on r_s is primarily through H(z), where the total radiation density Ω_r includes the neutrino contribution. The sound speed c_s = c/√(3(1+R_b)) depends only on the photon and baryon densities and does not require a neutrino correction (Hu Sugiyama 1995 [HuSugiyama1995]). The parameter Ω_Λ ≈ 0.699 is fixed by the bimetric model parameters; h and ω_b are the two free parameters in the CMB fit.
Applicability to modified gravity: The compressed CMB likelihood is derived under GR and the Planck Collaboration (2015, Results XIV) [Planck2015XIV] warns that it may be inapplicable to modified gravity models where perturbation growth differs from GR. We justify its use by noting that Vainshtein screening restores GR for perturbations at z_* = 1090, where H(z_*) ≫ m_FP and the bimetric correction is suppressed. Lüben et al. (2020) [Luben2020] used CMB data in their bimetric analysis, justifying it by the same Vainshtein screening argument. **Warning — perturbations not included.** The compressed likelihood probes the background geometry only. Beyond the shift parameters, bimetric gravity modifies the growth of structure and the decay of gravitational potentials throughout the unscreened region (z ≲ 12), sourcing an Integrated Sachs-Wolfe signal accumulated along the line of sight; a full Boltzmann treatment (CAMB/CLASS with the bimetric stress-energy tensor) could therefore constrain — or exclude — the model far more strongly than the background χ² suggests, and is the single most important follow-up. CMB E²(z): The full bimetric E²(z) is used for the D_M(z_*) integral. The ΛCDM E²(z) approximation, while adequate at θ = 6.1° (correction < 0.3%), underestimates Δχ²_CMB by up to 3.6× at intermediate θ and is not used.
### 3.4. Chi-Squared Calculation
For each probe, χ² = (data - model)^T C^(-1) (data - model), minimised over the relevant free parameters. For BAO, we marginalise over f = c/(H₀ r_d); for SNe, over M (and Ω_m for ΛCDM; the bimetric Ω_m is fixed); for CMB, over h and ω_b (and ω_Λ for ΛCDM; the bimetric Ω_Λ ≈ 0.699 is fixed by the model). The bimetric parameters (θ, m_FP, α, β) are held fixed at the baseline point θ = 6.1°, m_FP/H₀ = 25, α = 0, β = 3, which follows the best-fit configuration of Högås & Mörtsell (2025) for the DESI-era data combination; this choice anchors the comparison to the observationally favoured region of the parameter space. Radiation density: The BAO computation uses Ω_r = 5.38 × 10⁻⁵ (photons only, consistent with the DESI DR2 fiducial cosmology), while the SNe and CMB computations include the neutrino contribution (N_eff = 3.046, Ω_r = 9.08 × 10⁻⁵); this probe-dependent convention follows the fiducial treatment of each dataset, and the effect on Δχ² is < 0.01 for both probes, as Ω_r is negligible at z < 3. Numerical caveat: at the degenerate locus 1+2α+β = 0 (e.g., (α, β) = (−1, 1)) the leading coefficient of the quartic vanishes; for z ≳ 0.76 the solver then falls back to y = 1 (with a stderr warning), reverting E²(z) to its ΛCDM-like form. This point lies off the physical branch (B₁ = 0) and does not affect the canonical results, but grid values there should be read with this caveat.
## 4. Results
### 4.1. BAO
The bimetric model yields χ²_BAO = 10.464 vs. ΛCDM χ²_ΛCDM = 10.286, giving Δχ² = 0.178 (p = 0.66 for 1 degree of freedom, corresponding to the single effective parameter shift in Ω_m). We note that the Wilks theorem requires regularity conditions; at θ = 0, the parameter m_FP becomes unidentifiable, which technically violates these conditions. The p-values should therefore be interpreted as approximate. At the baseline parameters, BAO is effectively blind to all four bimetric parameters (i.e., insensitive at current data precision, and partially degenerate with the marginalised distance ratio f):
- (α, β) scan (20 combinations): Δχ² = 0.178 ± 0.001 --- flat across the entire grid.
- m_FP/H₀ scan (10 ≤ m_FP/H₀ ≤ 10⁵, four orders of magnitude): Δχ² ranges from 0.178 to 0.250, with the baseline value 0.178 at m_FP/H₀ = 25.
- θ scan: Δχ² grows from 0.178 at θ = 6.1° to ~2 at θ ≈ 31.5°, because the within-tracer DM-DH correlations in the block-diagonal covariance break the f-degeneracy at large angles. Our θ = 6.1° is deep in the blind regime (Δχ² ≈ 0.178), but BAO becomes weakly sensitive to θ at large angles. Under the ΛCDM null (10⁴ synthetic DESI DR2 realisations drawn from the full covariance, with f refit for both models), Δχ² ≥ 0.178 occurs with probability p = 0.52: the BAO penalty is statistically unremarkable.
The Δχ² decomposition reveals that the dominant contribution comes from the Ω_m shift (+0.662), which is partially compensated by the bimetric E²(z) correction (-0.485), yielding the net Δχ² = 0.178. The pure bimetric correction at fixed ΛCDM Ω_m contributes only Δχ² = 0.023; this is a standalone diagnostic and is not an additive term in the decomposition (0.662 - 0.485 = 0.177, while 0.662 - 0.485 + 0.023 = 0.200 ≠ 0.178). The bimetric E²(z) correction compensates for 73% of the Ω_m penalty (-0.485 out of +0.662). The 65× ratio of SNe pure correction (1.49) to BAO pure correction (0.023) arises because BAO observables are dimensionless ratios D_M/r_d and D_H/r_d, where the marginalisation over f = c/(H₀ r_d) absorbs the approximately constant component of the bimetric E²(z) correction, leaving only a residual z-dependent signal of order 0.023. SNe observables (luminosity distances D_L) have no analogous marginalisation: the absolute scale is set by H₀ through M, but the shape of D_L(z) directly reflects E²(z), so the full 1% correction contributes 1.49.
### 4.2. SNe Ia
Diagonal covariance. The bimetric model yields Δχ²_SNe = 15.58 for the baseline parameters, using the standard comparison: ΛCDM best-fit (Ω_m = 0.374, M marginalised) vs. bimetric (Ω_m = 0.304 from quartic, M marginalised), χ²/dof = 0.455. SNe Ia are not blind to the bimetric model:
- Standard comparison (diagonal): Δχ² = 15.58 (Ω_m^LCDM = 0.374 vs. Ω_m^bim = 0.304).
- Pure bimetric (same Ω_m = 0.304): Δχ² = 1.49---the E²(z) correction alone is sub-detectable.
- (α, β) scan: Δχ² = 13.9--15.6. m_FP/H₀ scan (≥ 5): Δχ² = 15.58. θ scan: Δχ² = 13.9--15.6 (flat, but non-zero).
Full covariance (STAT+SYS). Repeating the analysis with the full 1588 × 1588 Pantheon+ STAT+SYS covariance (analytical marginalization over M) gives: ΛCDM best-fit Ω_m = 0.347 (χ²_min = 1433.06, χ²/dof = 0.903); bimetric χ² = 1439.56; Δχ² = 6.50 (standard), and Δχ²_pure = 0.96. The full covariance thus reduces Δχ² by ~58%: the systematic correlations pull the ΛCDM best-fit Ω_m from 0.374 to 0.347---closer to the bimetric prediction (0.304) and to the BAO/CMB-preferred range---which shrinks the Ω_m mismatch that dominates Δχ². The decomposition mirrors the diagonal case: of the 6.50, 5.54 is attributable to the Ω_m shift and 0.96 to the pure E²(z) correction. The qualitative conclusion (SNe sensitive through the Ω_m prediction) is unchanged; the preference for ΛCDM is moderate: ~2.5σ with the full covariance, ~3.9σ with the diagonal approximation.
The Δχ² decomposition (diagonal case) reveals: Δχ² = -2 × cross-term + model-diff χ², where the cross-term measures the correlation of data residuals with the model difference, and model-diff χ² measures the centred shape difference. At Ω_m^LCDM = 0.374: model-diff = 15.64, cross-term = 0.03, giving Δχ² = 15.58. The Ω_m shift dominates: the bimetric E²(z) correction is negative (≲ 1%), which worsens the fit when Ω_m^bim < Ω_m^(LCDM best). For BAO, the same correction compensates (Ω_m^bim = 0.304 > Ω_m^(BAO best) = 0.297), but for SNe the direction is reversed. With the full covariance, this picture persists in weakened form: model shape differences contribute only 0.96, while the remaining 5.54 reflects the compressed Ω_m offset (0.347 - 0.304 = 0.043).
Intrinsic dispersion. To assess the robustness of the SNe sensitivity to the treatment of intrinsic scatter, we additionally compute Δχ²_SNe including a standard intrinsic dispersion σ_int = 0.12 mag, added in quadrature to the diagonal covariance (σ²_total = σ²_diag + σ²_int). This reduces Δχ²_SNe from 15.58 to 11.58 (Δχ²_pure from 1.49 to 1.11), while Ω_m^best remains 0.374. Even at σ_int = 0.20 mag (an upper bound), Δχ²_SNe = 8.06 --- SNe Ia remain a sensitive discriminator. The qualitative conclusion is unchanged. We do not add σ_int on top of the full STAT+SYS covariance: the intrinsic dispersion is already included in its diagonal (standard Pantheon+ construction), and adding it again would double-count the same physical effect. The σ_int scan (diagonal) is tabulated in Table A6.
The Ω_m = 0.374 obtained from the diagonal covariance is shifted 2.2σ from the Pantheon+ full-covariance value (0.334 ± 0.018 [Brout2022]). With full covariance, Ω_m^(LCDM best) = 0.347, closer to the bimetric prediction, reducing Δχ²; the Ω_m shift decreases from 0.070 to 0.043. This is a quantitative, not qualitative, change: SNe remain sensitive to the bimetric Ω_m prediction even with full covariance (6.50 > 3.84 at 95% CL).
### 4.3. CMB
The bimetric model yields χ²_CMB ≈ 0.04 (2 free parameters h, ω_b; Ω_Λ ≈ 0.699 is fixed by the model, 3 data points, DOF = 1), vs. ΛCDM χ² = 1.41 (2 free parameters h, ω_b with Ω_m and Ω_Λ fixed at the DESI-fiducial values; 3 data points, DOF = 1, p ≈ 0.24 — statistically consistent). The comparison thus rests on the relative constraining power (DOF = 1 for both models), not on absolute χ² saturation. We emphasise that Δχ²_CMB ≈ 0.04 should be interpreted as a diagnostic of model distance, not as a formal hypothesis test: both models are statistically consistent with the compressed data at baseline, and the discriminatory power comes from the θ-dependence, not from the baseline offset. For reference, Δχ² = 3.84 corresponds to the 95% confidence level for 1 DOF; thus even at the baseline θ = 6.1°, CMB mildly disfavours the bimetric model (Δχ² ≈ 0.02 ≪ 3.84 — a negligible offset), and does not exclude it. CMB is not blind to the mixing angle θ:
- Baseline (θ = 6.1°): Δχ² ≈ 0.04.
- (α, β) scan (20 combinations): Δχ² = 0.014--2.14. The predicted Ω_m is independent of (α, β) across the entire grid (0.3029 at all 20 points; |∂Ω_m/∂α| ≈ 4 × 10⁻⁹, |∂Ω_m/∂β| ≈ 10⁻¹², Appendix B.2), because the de Sitter attractor pins y(0) to 1 − 4.9 × 10⁻⁴ and the (α, β)-dependent terms in Ω_DE(0) enter only at O((1−y(0))²). Consequently, the SNe Ω_m tension cannot be alleviated by varying the potential parameters within this class: it is structural, not parametric. m_FP/H₀ scan (5 ≤ m_FP/H₀ ≤ 100): Δχ² = 0.001--0.36.
- θ scan: Δχ² grows from 0.04 to 2369 as θ increases from 6.1° to 30° (v10.6, with massive neutrinos evolved as matter at z ≲ 130; verified to <0.5% by an independent blind implementation). Note: at θ = 30° with m_FP/H₀ = 25, the dynamical Higuchi bound m_FP² > 2Ω_Λ (Högås & Mörtsell 2021, Eq. 5.11) is satisfied (625 > 1.4), so θ = 30° satisfies the static late-time bound for this non-self-accelerating model (the dynamical z-dependent Higuchi bound is stricter at high z and is noted as a caveat), though it is strongly excluded by the CMB geometric diagnostic (Δχ² = 2369).
The CMB-only sensitivity to θ, after isolating the BAO and SNe contributions, is presented for the first time in this work (Fig. 3). Högås & Mörtsell (2025) obtained θ constraints from the combined fit but did not decompose the constraining power by probe.
The θ grid (6.1°–30°) spans the posterior support of Högås & Mörtsell (2025) and the angles that would reproduce the externally preferred Ω_m values; the m_FP/H₀ scan samples four logarithmic decades. ### 4.4. Table 2: Summary of Δχ² scan results
```
Scan              BAO Δχ²      SNe Δχ²                CMB Δχ²
(alpha,beta), 20  0.178        13.9--15.6 (diag)      0.014--2.14
combos
m_FP/H0           0.178--0.250 15.58 (diag)           0.001--0.36
theta (baseline   0.178        15.58 (diag); 6.50     ~0.04
6.1deg)                        (full cov)
theta (up to      -> 2 at     13.9--15.6 (diag)       -> 2369
31.5deg)          31.5deg
Verdict           Blind        Sensitive               Sensitive to theta
```
Figure 3: Δχ² as a function of the mixing angle θ for each probe. BAO (blue) remains near zero at baseline θ = 6.1° but grows to ~2 at θ ≈ 31.5° as within-tracer correlations break the f-degeneracy, while SNe Ia (red) are offset at Δχ² ≈ 15 (diagonal; 6.5 with full covariance) and CMB (green) grows rapidly, reaching Δχ² 2369 at θ = 30°. The baseline θ = 6.1° is marked. The horizontal dotted line shows Δχ² = 3.84 (95% CL for 1 DOF).
## 5. Discussion
### 5.1. Vainshtein Screening: Per-Probe Mechanism
The blindness of BAO and the sensitivity of SNe and CMB arise from the interplay of two distinct mechanisms: the de Sitter attractor (y → 1) that suppresses the bimetric correction at late times, and the Vainshtein screening (Lüben et al. 2020 [Luben2020]) that operates in the early universe. The independence of these predictions across different redshift ranges confirms that the blindness and sensitivity are intrinsic dynamical properties of the model, not the result of parameter tuning.
The critical Hubble rate is H_c = m_FP (Lüben et al. 2020). For orientation, the Vainshtein radii implied by m_FP/H₀ = 25 are r_V(☉) ≈ 11 pc and r_V(Milky Way) ≈ 100 kpc; the GW170817 bound on the gravitational-wave speed is ~9 orders of magnitude weaker than would be needed to exclude this class, so screening neither rescues nor undermines the cosmological conclusions — it delineates their regime of applicability. For H > m_FP (early universe, z ≳ 12 at m_FP/H₀ = 25), Vainshtein screening is active and GR is restored. For H < m_FP (late universe, z ≲ 12), screening is inactive and the bimetric correction is present. In this regime, GR is approached not through Vainshtein screening but through the de Sitter attractor y → 1 (Lüben et al. 2020), which makes the correction small but non-zero.
BAO measures D_M/r_d and D_H/r_d at z < 2.33---the late universe where Vainshtein screening is inactive. The bimetric correction to E²(z) is present but intrinsically small (≲ 1.1% within the BAO range), because the dynamical evolution drives y → 1 at late times as a de Sitter attractor (y(0) = 0.99951), making Ω_DE ≈ Ω_Λ with only a tiny correction. The marginalisation over f = c/(H₀ r_d) further absorbs residual signal, because BAO observables are dimensionless ratios. The combination renders Δχ² ≈ 0.
SNe Ia measure D_L(z) = (1+z) D_M(z) at z < 2.26---the same late-time regime where screening is inactive. The H(z) integral accumulates only the small unscreened correction (≲ 1.2%), and the pure bimetric correction at fixed Ω_m gives Δχ² = 1.49 (diagonal) / 0.96 (full covariance)---sub-detectable. However, SNe Ia are sensitive to the bimetric model through the Ω_m shift: the quartic constraint predicts Ω_m = 0.304, while ΛCDM best-fit to Pantheon+ gives Ω_m = 0.374 (diagonal) / 0.347 (full). The resulting Δχ² = 15.58 / 6.50 is dominated by this Ω_m mismatch, not by the E²(z) correction itself. Unlike BAO, there is no r_d marginalisation to absorb the signal: D_L depends directly on H₀ and Ω_m through the luminosity distance integral. The negative bimetric E²(z) correction shifts the model in the wrong direction relative to the SNe-preferred Ω_m, amplifying rather than compensating the Ω_m discrepancy. This contrasts with BAO, where Ω_m^bim = 0.304 > Ω_m^(BAO best) = 0.297 and the correction partially compensates.
CMB measures D_M(z_*)---the comoving distance integrated from z = 0 to z_* = 1090---and the sound horizon r_s(z_*) = ∫_{z*}^∞ c_s(z')/H(z') dz' (from z_* to ∞). The integral D_M(z_*) traverses two regimes: the unscreened late universe (z ≲ 12, where H < m_FP and the bimetric correction is present) and the screened early universe (z ≳ 12, where H > m_FP and the correction is suppressed). Although the correction in the unscreened regime is small (≲ 1.2%), the unscreened region (z ≲ 12) contributes 72% of D_M(z_*) (Table A4; computed as ∫₀¹² dz/E(z) / ∫₀¹⁰⁹⁰ dz/E(z) ≈ 0.72 for the bimetric model with baseline parameters), so the correction accumulates and produces a detectable shift in D_M(z_*), which the compressed CMB likelihood registers as Δχ² > 0. The sound horizon r_s(z_*), being an integral from z_* to ∞, lies entirely in the screened regime and is less affected. This asymmetry between the screened early universe and the unscreened late universe is the fundamental reason CMB is a sensitive discriminator among the three probes, alongside SNe which are sensitive through the Ω_m prediction.
**Table 3: Per-probe blindness and sensitivity mechanism summary**
```
Probe    Bimetric        Marginalisation     Source of sensitivity/blindness
         correction to                          
         E2(z)                                 
BAO      ~< 1.1% (z<=2.33) Yes: f = c/(H0 rd)   Blind: small correction +
                                             f-marginalisation absorbs
                                             residual
SNe Ia   ~< 1.1% (z<=2.26) No                 Sensitive: Omega_m mismatch
                                             (0.304 vs 0.374 diag / 0.347
                                             full) dominates (~90% diag,
                                             ~85% full)
CMB (D_M) ~< 1.2% (z ~< 12, No               Sensitive: 72% of D_M(z*) from
          unscreened)                          unscreened region accumulates
                                             correction
CMB (r_s) Suppressed (z >=  No                Blind: r_s integral entirely in
          1090, screened)                      screened regime
```
### 5.2. Phantom Dark Energy Phase
The effective equation of state w_DE(z) < -1 (phantom) at all redshifts where Ω_DE is significant (z ≲ 5) is a genuine prediction of the bimetric model. BAO is insensitive to this phase because the bimetric correction to E²(z) remains small (≲ 1.1% within the BAO redshift range (z ≤ 2.33), rising to ~1.09% at z ≈ 3 and ~1.12% at z ≈ 5 just beyond the BAO range, for context). SNe are sensitive to the bimetric model not through the E²(z) correction directly (Δχ²_pure = 1.49 / 0.96), but through the Ω_m prediction (Δχ² = 15.58 / 6.50). The CPL approximation (w_0 ≈ -0.861, w_a ≈ -0.500) fitted on z [0, 3] does not capture the phantom nature (w_DE(0) = -1.005); Ω_DE(z) crosses zero at z ≈ 5.1 and tends to the negative asymptote Ω_Λ − s²m²(1+β/3) = −13.42 (Ω_DE(10) ≈ −3.11), making CPL inapplicable beyond z ≈ 4. The discrepancy arises because the effective w_DE(z) in the bimetric model has a complex redshift dependence---with Ω_DE(z) changing sign at z ≈ 5.1 and becoming negative at higher z---that a 2-parameter CPL fit cannot capture.
### 5.3. Combined Δχ² and the Dominance of SNe
If the three probes are combined, the total Δχ² is dominated by SNe. To keep the combination statistically interpretable, we quote two internally consistent schemes rather than mixing covariance choices across probes:
```
Scheme          BAO      SNe      CMB      Sum Δχ²
All-diagonal    0.195    15.58    0.45     16.22
All-full        0.178    6.50     0.04     6.72
(also: σ_int = 0.12, diagonal SNe scheme: 0.195 + 11.58 + 0.45 = 12.22)
```
Here "all-diagonal" uses the block-diagonal BAO result with rounded Table-IV uncertainties (0.195), the diagonal SNe covariance (15.58), and the CMB value from the LCDM-E²-based recalculation (0.45); "all-full" uses the full 13×13 BAO covariance (0.178), the full Pantheon+ STAT+SYS SNe covariance (6.50), and the bimetric-E² CMB θ-scan at baseline (0.04, v10.6 with massive-neutrino evolution). A hybrid combination (full-covariance BAO with diagonal SNe) mixes error models across probes and is not statistically interpretable; we do not quote it. We assume the three probes are statistically independent, which is a good approximation for BAO, SNe Ia, and CMB distance priors; residual correlations could arise through shared nuisance parameters (e.g., the absolute distance scale), and a full joint analysis with a block-diagonal super-covariance is left to the MCMC extension of this work.
SNe contribute ~96% of the total at baseline in the all-diagonal scheme and ~97% in the all-full scheme, through the Ω_m mismatch. CMB contributes ~0.3% at baseline but grows rapidly with θ. In the combined analysis of Högås & Mörtsell (2025), both SNe and CMB drive the constraints on θ, while BAO contributes little.
### 5.4. Implications for Future Surveys
The bimetric correction to E²(z) at BAO redshifts is small (≲ 1.1%) but non-zero---Vainshtein screening is inactive at z < 12, so the correction is present, not suppressed. The current BAO blindness therefore reflects the limited precision of existing data (~1--2%) combined with the marginalisation over f, not a fundamental theoretical limitation. Future sub-percent BAO surveys (DESI DR3, Euclid) could in principle detect the z-dependent component of the 1% correction, though the f-marginalisation absorbs any approximately constant (in z) component, leaving only the residual z-dependent signal as potentially detectable. Probes that integrate over longer paths, such as CMB, or which directly sample the transition region — 21 cm cosmology, including the global 21-cm signal at z ≈ 15–30 — remain the most promising avenue for constraining θ.
### 5.5. Ω_m and H₀
We note that the predicted Ω_m = 0.304 is in mild tension with independent measurements: Planck 2018 gives Ω_m = 0.315, and the Pantheon+ full-covariance best-fit gives Ω_m = 0.334 ± 0.018. This tension is itself a constraint on the bimetric model at the baseline parameter point. Two structural results sharpen it into a no-go statement. First, the (α, β) scans (§4.3, Table A1b) show that the predicted Ω_m is independent of the potential parameters to |∂Ω_m/∂α| ≈ 4 × 10⁻⁹, |∂Ω_m/∂β| ≈ 10⁻¹²: within this class of interactions, no choice of (α, β) can move Ω_m toward the data-preferred 0.347. The only remaining dials are (θ, m_FP, Ω_Λ), and each is CMB-taxed: raising θ to reach Ω_m = 0.315 costs Δχ²_CMB ≈ 32. Varying the anchor Ω_Λ fares no better: along the degeneracy ray (closure Ω_m·cos²θ = 1 − Ω_Λ − Ω_r), at Ω_Λ = 0.657 the model reaches the SNe best-fit region (Ω_m = 0.347, Δχ²_SNe = 0.05), but the CMB Δχ² rises by ≈ 28.6 and the BAO Δχ² by ≈ 23.4 (Δχ²_BAO(Ω_Λ) has a sharp minimum at the baseline anchor: +23.6 at Ω_Λ = 0.657, +6.3 at 0.680, +4.2 at 0.720, +21.4 at 0.740). The combined cost of moving along the ray is ≈ 45.6: the anchor is a sharp joint minimum of the three probes, not a hiding spot. Second, the SNe tension persists at Δχ² ≈ 13–16 across the entire (α, β) grid. The Ω_m tension is therefore structural, not parametric: either the baseline anchor (θ, m_FP, Ω_Λ) is excluded, or the Ω_m difference is a falsifiable prediction of the model. A full exploration of the (θ, m_FP, Ω_Λ) space is left for future work.
The self-consistent Ω_m = 0.304 differs from both the Planck value (0.315) and the DESI BAO best-fit (0.2972). This is not a discrepancy but a prediction: the quartic constraint determines Ω_m from the bimetric parameters, not from data. The H₀ = 68.8 km/s/Mpc follows from the BAO distance fit with f = 29.607 and Planck r_d = 147.09 Mpc (Planck 2018 TT,TE,EE+lowE best-fit [Planck2018], consistent with the DESI DR2 fiducial value r_d = 147.05 Mpc adjusted for ω_b = 0.02237), consistent with the joint analysis of Högås & Mörtsell (2025), who obtained H₀ = 69.0 ± 0.4 km/s/Mpc. The bimetric model does not claim to resolve the Hubble tension in this analysis: Lüben et al. (2020) [Luben2020] showed that resolution requires combining CMB, SNe Ia, and BAO data, where the phantom phase in the z ~1--10 interval shifts the inferred H₀ toward the local value. We note that Högås & Mörtsell (2026) [Hogas2026] have shown that physically motivated priors in the local distance ladder reduce the Hubble tension from 5σ to 2σ within ΛCDM.
To quantify whether the Ω_m tension can be alleviated within the bimetric parameter space, we scan Ω_m as a function of the mixing angle θ (fixing α = 0, β = 3, m_FP/H₀ = 25). Increasing θ from 6.1° to 20° raises the predicted Ω_m from 0.304 to 0.341: the larger mixing angle shifts y(0) further from unity, which increases the subtracted massive-graviton term and thereby decreases the effective dark-energy density, raising Ω_m = 1 - Ω_r - Ω_DE(0). However, this alleviation comes at a steep cost: the CMB Δχ² grows from 0.04 to 370 over the same range (Table A5). The Planck value Ω_m = 0.315 is reached at θ ≈ 12.2°, where Δχ²_CMB ≈ 31.6 (> 99% CL), and the Pantheon+ full-covariance value Ω_m = 0.334 requires θ ≈ 18.3°, where Δχ²_CMB ≈ 245. This indicates that the Ω_m tension is not a mild statistical offset but a structural feature of the bimetric model: reconciling Ω_m with external estimates requires mixing angles that are strongly excluded by CMB data.
### 5.6. BIC Comparison
The Bayesian Information Criterion for BAO gives ΔBIC_min = 0.178 (from the (α, β) scan minimum, matching but not derived from the baseline Δχ² = 0.178) (bimetric parameters fixed, δ k = 0)---models are indistinguishable for BAO. For the full model (δ k = 4), ΔBIC_full = 0.178 + 4ln(13) = 10.4, favouring ΛCDM. This BIC applies to BAO only. In the BAO analysis, all four bimetric parameters are held fixed; the BIC with δk = 4 is a hypothetical estimate for the full model, not a fit to BAO data. For SNe, Δχ² = 15.58 (6.50 with the full covariance) disfavours the bimetric model at the baseline point, though the diagonal-covariance value is partly driven by the diagonal covariance approximation. The CMB contribution (Δχ² ≈ 0.04) is negligible at baseline (Δχ² < 3.84 at 95% CL for θ = 6.1°).
## 6. Limitations and Future Work
This analysis has several limitations:
1. Diagonal SNe covariance (historical upper bound): The baseline Δχ²_SNe = 15.58 uses the diagonal approximation. With the full Pantheon+ STAT+SYS covariance (§4.2), Δχ²_SNe = 6.50 with Ω_m^(LCDM best) = 0.347: the tension is reduced but not eliminated. Both estimates are reported; 15.58 is an upper bound.
2. SNe intrinsic dispersion: The baseline diagonal analysis does not include an intrinsic dispersion σ_int in the SNe covariance. Including σ_int = 0.12 mag reduces Δχ²_SNe from 15.58 to 11.58 but does not change the qualitative conclusion (see §4.2, Appendix A Table A6). The baseline value Δχ²_SNe = 15.58 (diagonal) should be interpreted as an upper bound on the SNe sensitivity. σ_int is not added to the full STAT+SYS covariance (double counting).
3. Compressed CMB likelihood: The Bansal & Huterer (2025) compressed likelihood is derived under GR. Its applicability to bimetric gravity relies on Vainshtein screening restoring GR at z_* = 1090. This assumption should be tested with full Boltzmann codes (CAMB/CLASS) modified for bimetric gravity. We estimate that perturbative corrections to the compressed likelihood at z_* = 1090 (where Vainshtein screening suppresses the bimetric modification by ~ (H₀/m_FP)² ~ 0.16%) would shift Δχ²_CMB by ≲ 10^-3, well below the baseline value of 0.50, though a definitive assessment requires a full Boltzmann computation.
4. An updated compressed CMB likelihood, CMBComp [Srivastav2026] (Srivastav, Bansal & Huterer 2026, arXiv:2606.18455), incorporating SPT-3G D1, ACT DR6, Planck PR3 primary CMB anisotropy, and Planck PR4/NPIPE CMB-lensing data and specifically designed for modified expansion histories, is now available; adopting it is a natural extension for future work.
5. Grid scan over θ: a full MCMC analysis over (θ, m_FP, α, β) would provide proper posterior distributions and covariances. The θ-scan Δχ² values should be interpreted as a sensitivity map rather than independent hypothesis tests; a trial-corrected (look-elsewhere) significance would be needed for a formal detection claim at a specific θ.
6. Fixed α and β: only α = 0, β = 3 are tested at baseline. The (α, β) scans show Δχ² is flat for BAO (near zero) and flat but non-zero for SNe (14--16), while CMB shows mild variation (~0.50--0.6). The untested parameter space includes: α ≠ 0 (only 20 grid points scanned), β ≠ 3, combined (θ, m_FP) variations (scanned independently, not jointly), and m_FP/H₀ < 10 (not explored). Conclusions about BAO blindness apply to the tested grid, not to the full parameter space.
7. Neutrino treatment: massive neutrinos (Σm_ν = 0.06 eV) are evolved as nonrelativistic matter at z ≲ 130 and as radiation above (Bansal & Huterer 2025, Appendix C convention), and the R prefactor follows their Eq. 4; the sound speed c_s = c/√(3(1+R_b)) carries no neutrino correction, following Hu Sugiyama (1995). The residual convention uncertainty was tested numerically: varying the neutrino transition redshift over Z_NR ∈ [100, 200] shifts Δχ²_CMB(6.1°) by only 0.018--0.030, and the ν-matter treatment shifts intermediate-θ values by up to ~15% (largest at small χ²); conclusions are unchanged.
8. DESI DR2 only: future data releases may tighten BAO constraints; the bimetric correction is present at z < 12 (screening inactive), so improved precision could in principle detect it.
14. Blindness is not an artefact of f-marginalisation: the SNe probe involves no r_d-marginalisation, yet its pure bimetric Δχ² is 0.96 — the de Sitter attractor (ε = 4.9 × 10⁻⁴) suppresses the correction before any marginalisation, and f absorbs only the constant r_d rescaling, not z-dependent distortions.
15. A full joint MCMC over (θ, m_FP, Ω_Λ) with an emulator remains future work; the present per-probe decomposition is exact at the anchor and along the degeneracy ray (§5.5), but does not replace posterior exploration.
Future work should:
9. Test the compressed CMB likelihood with full Boltzmann codes for bimetric gravity.
10. [Addressed in this work, §4.2: the full Pantheon+ covariance gives Δχ²_SNe = 6.50.] Remaining extension: a full MCMC over (θ, m_FP, α, β) for the combined BAO + SNe + CMB likelihood.
11. Investigate whether m_FP/H₀ < 10 produces a threshold where BAO becomes sensitive.
12. Explore whether varying α and β can break the flatness of Δχ² for BAO.
13. Test predictions with 21 cm cosmology, which may probe the transition between the screened early universe and the unscreened late universe (z ≲ 12).
## 7. Conclusions
We have performed a probe-by-probe test of bimetric gravity using DESI DR2 BAO, Pantheon+ SNe Ia, and compressed Planck 2018 + ACT DR6 data. Our main results are:
15. BAO is effectively blind to bimetric parameters at baseline: Δχ² = 0.178 (p ≈ 0.66, 1 DOF; Wilks conditions violated at θ = 0, see §4.1), flat across 20 (α, β) combinations (Δχ² = 0.178 ± 0.001) and varying only from 0.178 to 0.250 across four orders of magnitude in m_FP/H₀ (10 ≤ m_FP/H₀ ≤ 10⁵); BAO becomes sensitive to θ only at large angles (Δχ² → 2 at θ ≈ 31.5°, see §4.1).
16. SNe Ia are not blind to the bimetric model: Δχ² = 15.58 (diagonal covariance --- an upper bound), 11.58 with σ_int = 0.12 mag, and 6.50 with the full Pantheon+ STAT+SYS covariance --- a moderate, 2.5σ--3.9σ preference for ΛCDM depending on the error model. With full covariance the ΛCDM best-fit Ω_m shifts from 0.374 to 0.347, closer to the bimetric prediction (0.304) and to the BAO/CMB-preferred range, which accounts for most of the reduction (the Ω_m mismatch remains the dominant contribution: 5.54 of 6.50, vs 0.96 for the pure E²(z) correction). The predicted Ω_m is independent of the potential parameters (α, β) to ≲ 10⁻⁸ (Table A1b), so this tension is structural and cannot be relieved by reshaping the interaction potential.
17. CMB is a sensitive discriminator: Δχ² ≈ 0.04 at baseline, growing to 2369 as θ → 30°. CMB and SNe both constrain the bimetric model; BAO does not.
18. The BAO blindness mechanism is the smallness of the bimetric correction at late times: BAO probes only z ≲ 2.3, where the de Sitter attractor (y → 1) renders the correction ≲ 1.1%, and marginalisation over f absorbs residual signal. SNe probe the same redshift range but lack an analogous marginalisation: the Ω_m prediction from the quartic constraint is directly tested, yielding Δχ² = 15.58 (diagonal; upper bound) / 6.50 (full). Vainshtein screening is inactive in this regime (H < m_FP). CMB integrates D_M(z_*) from z = 0 to z_* = 1090, accumulating the unscreened late-universe correction (z ≲ 12, contributing 72% of D_M(z_*)), while the screened early universe (z ≳ 12) contributes little.
19. The pure bimetric correction at fixed ΛCDM Ω_m contributes only Δχ² = 0.023 for BAO and Δχ² = 1.49 (diagonal) / 0.96 (full covariance) for SNe; the net BAO Δχ² = 0.178 is dominated by the Ω_m shift, while the net SNe Δχ² (15.58 / 6.50) is also dominated by the Ω_m shift, but in the opposite direction.
20. The bimetric model predicts a phantom dark energy component (w_DE < -1) at all redshifts where Ω_DE is significant (z ≲ 5), with w_DE(0) = -1.005.
21. The bimetric correction at BAO redshifts is small (≲ 1.1%) but present (Vainshtein screening is inactive at z < 12). Current BAO precision (~1--2%) is insufficient to detect it, but future sub-percent surveys could in principle detect the z-dependent residual, subject to the f-degeneracy absorbing any approximately constant component.
22. BIC: ΔBIC_min = 0.178 (from the (α, β) scan minimum, matching but not derived from the baseline Δχ² = 0.178) (indistinguishable at fixed parameters), ΔBIC_full = 10.4 (ΛCDM preferred with 4 free bimetric parameters).
This work reproduces and extends the results of Högås & Mörtsell (2025) [Hogas2025] by providing an independent per-probe pipeline, explicit B-parameter values, the per-probe Δχ² decomposition, the CMB-only θ sensitivity, and the full-covariance SNe assessment (Δχ² = 6.50). The key finding is that BAO is blind to bimetric gravity (Δχ² = 0.178), while SNe Ia (Δχ² = 6.50--15.58 depending on the covariance) and CMB (Δχ² ≈ 0.04 at baseline, growing steeply with θ) are both sensitive. The phantom effective dark energy (w_DE < -1) of the model qualitatively moves in the direction of easing both the H₀ and S₈ tensions, as found by Lüben et al. (2020) in the joint analysis; a quantitative assessment requires the Boltzmann-level extension. Future work should combine all three probes with full covariance matrices and MCMC, and test the compressed CMB likelihood with bimetric-modified Boltzmann codes.
## Acknowledgments
The computational work in this study was performed with the assistance of AI tools. The early main numerical computations were carried out with Alice AI (Yandex). The current phase — administration of the project workflow, code development, and numerical computations, including the quartic solver, χ² minimisation, B-parameter evaluation, and the full-covariance supernova pipeline — was carried out using Python code developed in collaboration with Kimi (Moonshot AI), in the role of the project's AI administrator ("Алиса"). The manuscript underwent multiple rounds of review by AI-assisted auditors (GigaChat Standard, GigaChat Reasoning, and GigaChat Research, Sber), whose feedback identified formatting issues, logical inconsistencies, and bibliographic errors. Independent computational cross-checks of key results were performed by additional AI systems used as verifiers: DeepSeek (DeepSeek AI) confirmed all seven verified quantities (B-parameters, quartic root, Ω_m self-consistency, H₀, ΔBIC, and Δχ² decompositions) to full precision, and — decisively for the CMB results — the compressed-likelihood θ-scan was independently reproduced by three separate, self-contained codebases written by DeepSeek, Qwen (Alibaba), and GLM (Zhipu AI), each executed verbatim by the project administrator; all three reproduced the canonical values within 0.1% (with one of these cross-checks localising a real sign error in the sound-speed formula, fixed in archive v10.6). The full-covariance supernova analysis was independently reproduced from the public Pantheon+ data before inclusion (§4.2). All AI-generated content was independently verified by the author; in particular, numerical claims supplied by external AI systems were accepted only after verbatim execution of their code, never on the basis of reported values alone. The author takes full responsibility for all results and conclusions presented in this work.
## Data Availability
The DESI DR2 BAO data are publicly available at https://data.desi.lbl.gov/doc/releases/dr2/ [DESI2025a, DESI2025b]. The Pantheon+ data are publicly available at https://github.com/PantheonPlusSH0ES/DataRelease [Brout2022, Scolnic2022]. The y(z) table, B-parameter values, w_DE(z) values, and Δχ² scan results are provided in this paper. The Python code (archive v10.6) is publicly available at https://github.com/chernicovalex-prog/bao-blindness-bimetric and archived with DOI 10.5281/zenodo.23292160.
## References
[Hassan2012] S. F. Hassan and R. A. Rosen, "Bimetric Gravity from Ghost-Free Massive Gravity," JHEP 02, 126 (2012), arXiv:1109.3515.
[vonStrauss2012] M. von Strauss, A. Schmidt-May, J. Enander, E. Mörtsell, and S. F. Hassan, "Cosmological Solutions in Bimetric Gravity and their Observational Tests," JCAP 03, 042 (2012), arXiv:1111.1655.
[Comelli2012] D. Comelli, M. Crisostomi, F. Nesti, and L. Pilo, "FRW Cosmology in Ghost Free Massive Gravity," JHEP 03, 067 (2012), arXiv:1111.1983.
[Hogas2021] F. Högås and E. Mörtsell, "Analytical constraints on bimetric gravity," JCAP 05, 001 (2021), arXiv:2101.08794.
[Hogas2021b] F. Högås and E. Mörtsell, "Observational Constraints on Bimetric Gravity," JCAP 05, 002 (2021), arXiv:2101.08795.
[Hogas2025] F. Högås and E. Mörtsell, "Bimetric gravity improves the fit to DESI BAO and eases the Hubble tension," Phys. Rev. D 112, 103515 (2025), arXiv:2507.03743.
[Luben2020] M. Lüben, A. Schmidt-May, and J. Smirnov, "Vainshtein Screening in Bimetric Cosmology," Phys. Rev. D 102, 123529 (2020), arXiv:1912.09449.
[DESI2025a] DESI Collaboration, "DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints," Phys. Rev. D 112, 083515 (2025), arXiv:2503.14738.
[DESI2025b] DESI Collaboration, "DESI DR2 Results I: Baryon Acoustic Oscillations from the Lyman Alpha Forest," Phys. Rev. D 112, 083514 (2025), arXiv:2503.14739.
[Eisenstein2005] D. J. Eisenstein et al., "Detection of the Baryon Acoustic Peak in the Large-Scale Correlation Function of SDSS Luminous Red Galaxies," Astrophys. J. 633, 560 (2005), arXiv:astro-ph/0501171.
[Brout2022] D. Brout et al., "The Pantheon+ Analysis: Cosmological Constraints," Astrophys. J. 938, 110 (2022), arXiv:2202.04077.
[Scolnic2022] D. Scolnic et al., "The Pantheon+ Analysis: The Full Data Set and Light-curve Release," Astrophys. J. 938, 113 (2022), arXiv:2112.03863.
[Bansal2025] P. Bansal and D. Huterer, "Expansion-history preferences of DESI DR2 and external data," Phys. Rev. D 112, 023528 (2025), arXiv:2502.07185.
[HuSugiyama1995] W. Hu and N. Sugiyama, "Small Scale Cosmological Perturbations: An Analytic Approach," Astrophys. J. 444, 489 (1995), arXiv:astro-ph/9407093.
[Hogas2026] F. Högås and E. Mörtsell, "Physically-motivated priors in the local distance ladder significantly reduce the Hubble tension," MNRAS 548, stag724 (2026), arXiv:2601.22215.
[Planck2015XIV] Planck Collaboration, "Planck 2015 results. XIV. Dark energy and modified gravity," Astron. Astrophys. 594, A14 (2016), arXiv:1502.01590.
[Planck2018] Planck Collaboration, "Planck 2018 results. VI. Cosmological Parameters," Astron. Astrophys. 641, A6 (2020), arXiv:1807.06209.
[Srivastav2026] S. Srivastav, P. Bansal, and D. Huterer, "CMBComp: A compressed CMB likelihood for modified expansion histories," arXiv:2606.18455 (2026).
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

### Table A2: y(z) and E²(z) --- Bimetric vs ΛCDM
Fractional correction δE² = (E²_bim − E²_LCDM)/E²_LCDM. Maximum |δE²| ≈ 1.12% at z ≈ 5. The correction is negative (bimetric expansion is slower) at all z > 0.
```
z       y(z)      E2_bimetric    E2_LCDM     dE2 (%)
0.00    0.999512  1.000000       1.000000    +0.0000
0.10    0.999350  1.099643       1.100795    -0.1047
0.50    0.998352  1.715027       1.723295    -0.4798
1.00    0.996093  3.107718       3.132092    -0.7782
2.00    0.986814  8.830855       8.921404    -1.0150
2.33    0.981968  11.821478      11.946601   -1.0473
3.00    0.968770  19.980325      20.199724   -1.0862
5.00    0.896071  65.815888      66.561436   -1.1201
10.00   0.556527  403.379176     407.168012  -0.9305
12.00   0.442652  666.804264     672.033688  -0.7781
1090    0.000002  523876594.328  523922423.336 -0.0087
```
### Table A3: Effective Dark Energy Equation of State w_DE(z)
Phantom regime (w_DE < −1) at all z where Ω_DE is significant (z ≲ 5). Ω_DE → 0 at z ≈ 5, making w_DE ill-defined beyond that point.
```
z       w_DE(z)
0.0     -1.0050
0.5     -1.0169
1.0     -1.0411
2.0     -1.1536
3.0     -1.4610
4.0     -2.5988
4.5     -5.4906   (ΩDE -> 0; ill-defined)
5.0     16.2874   (ΩDE -> 0; ill-defined)
```
### Table A4: D_M(z*) Accumulation by Redshift
Fraction of the total comoving distance D_M(z* = 1090) contributed by z ≤ z_check. The unscreened region (z ≲ 12, where H < m_FP) contributes 72.19% of D_M(z*), despite spanning only ~1% of the redshift range to z*.
```
z <= z_check    % of D_M(z*)    Regime
1               24.37%          Unscreened
5               57.20%          Unscreened
10              69.43%          Unscreened
12              72.19%          Screening boundary
50              87.85%          Screened
500             98.51%          Screened
1090            100.00%         z* (last scattering)
```
### Table A5: CMB θ-Scan --- Δχ² vs Mixing Angle
CMB-only Δχ² (BAO and SNe contributions isolated), computed with the v10.6 pipeline (massive neutrinos as matter at z ≲ 130), verified to <0.5% by an independent blind implementation (§4.3). At θ = 30° with m_FP/H₀ = 25, the Higuchi bound is satisfied (625 > 1.4).
```
θ (deg)    Δχ²_CMB    Significance
6.1        0.04       < 95% CL
10.0       9.3        > 95% CL
12.2       31.6       > 99% CL
15.0        94        >> 99% CL
18.3       245        >> 99% CL
20.0       370        >> 99% CL
25.0       1029       >> 99% CL
30.0       2369       >> 99% CL
```
### Table A6: SNe σ_int Scan --- Δχ² vs Intrinsic Dispersion (diagonal covariance)
Diagonal covariance with σ_int added in quadrature: σ²_total = σ²_diag + σ²_int. Ω_m^best refers to the ΛCDM fit. Even at σ_int = 0.20 mag (upper bound), Δχ²_SNe = 8.06 --- SNe remain a sensitive discriminator.
```
σ_int (mag)   Ω_m^best (ΛCDM)   Δχ² (standard)   Δχ² (pure)   χ²/dof
0.00          0.3742             15.58            1.49         0.455
0.05          0.3744             14.68            1.40         0.428
0.08          0.3744             13.49            1.29         0.393
0.10          0.3744             12.55            1.20         0.365
0.12          0.3742             11.58            1.11         0.338
0.15          0.3740             10.15            0.98         0.298
0.20          0.3736             8.06             0.78         0.239
```
## Appendix B.2: Analytic Late-Time Expansion and the Rigidity of Ω_m

At z = 0 the de Sitter attractor pins y to unity: writing y = 1 − ε, the quartic
F(ε) ≡ f(1−ε; Ω_m(ε)) with the closure Ω_m(ε) = 1 − Ω_r − ω_ν/h² − Ω_Λ +
sin²θ·m²_FP·ε·(1 + αε + βε²/3) expanded to first order gives

  F(0) = 1 − Ω_Λ − Ω_r − ω_ν/h²          (≈ 0.2996 at baseline)
  F'(0) = −[f'(1) − sin²θ·m²_FP]          (≈ −617 at baseline, dominated by the
                                           spin-2 mass scale m²_FP = 625)
  ε = F(0)/|F'(0)| + O(ε²)                (= 4.857 × 10⁻⁴ vs the exact
                                           4.859 × 10⁻⁴ — the linear expansion
                                           reproduces ε to 0.05%)

so ε ∝ (1 − Ω_Λ − Ω_r − ω_ν/h²)/m²_FP is set by the distance from the de Sitter
point, not by the potential shape. The (α, β)-dependence of Ω_DE(0) is
explicitly O(α·ε²) and O(β·ε³), and the implicit dependence through ε adds only
|∂Ω_m/∂α| ≈ 4 × 10⁻⁹ and |∂Ω_m/∂β| ≈ 10⁻¹² (finite differences of the exact
quartic solution at baseline). Hence the present-day matter density is a
function of (θ, m_FP, Ω_Λ) alone within this model class. Note that the rigidity
is an amplitude statement at z = 0: at z ≳ 0.5 the deviation 1−y(z) grows and
(α, β) do shape E²(z) — they are constrained through the form of the expansion
history (the pure corrections of §4), not through the present-day energy budget.

## Appendix B: Pseudocode for E²(0) Residual Computation
The E²(0) residual (~5 × 10⁻¹⁵) serves as a numerical consistency check. The following pseudocode describes the computation:
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
     c. Solve a₀ + a₁y + a₂y² + a₃y³ + a₄y⁴ = 0 using Brent's method on [0, 1.0]
     d. Compute E²(0) = Ω_m + Ω_r + Ω_DE(y(0)),
        Ω_DE = Ω_Λ − sin²θ·m²_FP·(1−y)·[1 + α(1−y) + β/3·(1−y)²]
   Find Ω_m such that E²(0) = 1 (Brent's method)
3. Compute residual: E²(0)_residual = |E²(0) − 1.0|; expected ~5e-15
4. Round-trip test: (θ, m_FP, α, β, Ω_Λ) -> B -> y -> recover params; ~1e-13
```
*--- End of Appendix ---*
