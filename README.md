# BAO Blindness to Bimetric Gravity — Reproducible Code

Companion code to the paper *"BAO Blindness to Bimetric Gravity:
A Per-Probe Δχ² Decomposition"* (A. Chernikov, 2026, v21.5).

The repository reproduces every canonical number of the paper from scratch.
Physics: Hassan–Rosen ghost-free bimetric gravity in the Högås & Mörtsell
B-parameter formalism; baseline θ = 6.1°, m_FP/H₀ = 25, α = 0, β = 3,
Ω_Λ = 0.699 (Högås & Mörtsell 2025, arXiv:2507.03743).

## Quick start (5 minutes)

```bash
pip install -r requirements.txt
make test        # self-checks + regression vs canonical numbers (~1-2 min)
make bao         # BAO: DESI DR2, Δχ² = 0.178 (full covariance)
make cmb         # CMB θ-scan, 8 angles (several minutes)
bash data/download_data.sh && make sne-full   # SNe full covariance (6.50)
```

## Layout

| Path | Content |
|---|---|
| `bimetric/core.py` | All background physics: quartic for y = a_f/a_g (H&M 2025, A13), E²(z) with the massive-ν transition at 1+z = 131, sound horizon, D_M, B-parameters. One implementation used by every probe. |
| `probes/bao.py` | DESI DR2 (13 points), full 13×13 covariance or Table IV block-diagonal |
| `probes/cmb_theta_scan.py` | Compressed Planck 2018 + ACT DR6 likelihood (Bansal & Huterer 2025), θ-scan 6.1°–30° |
| `probes/sne.py` | Pantheon+ (1588 SNe, z > 0.01), diagonal or full STAT+SYS covariance, analytic M-marginalisation |
| `canonical_numbers.json` | Every canonical number + tolerance + verification status |
| `tests/run_checks.py` | Self-checks (quartic identity, Brent vs `numpy.roots`, …) + regression against the canonical numbers |
| `verification/` | Independent verification artefacts: third-party code (written by external AI systems), run verbatim by the project, with logs and checksums |

## Canonical results (archive v10.6)

| Probe | Quantity | Value |
|---|---|---|
| Baseline | y(0), Ω_m, E²(0)−1 | 0.999512, 0.30439, 5·10⁻¹⁵ |
| BAO | Δχ² (full cov) | **0.178** (p ≈ 0.66 — BAO is blind) |
| SNe | Δχ² diagonal / full | **15.58** (upper bound) / **6.50** |
| CMB | Δχ² at θ = 6.1° / 30° | **0.04** / **2368.93** |
| Combined | all-full / all-diag | **6.72** / **16.22** |

Statistical conventions: σ_int is never added on top of the full STAT+SYS
covariance (double counting); the diagonal SNe value is an upper bound.

## Verification history

The CMB θ-scan was reproduced by **three independent codebases**, written by
external AI systems from a blind specification (no target values) and executed
verbatim by the project administrator:

| Code | χ²(6.1°) | χ²(15°) | χ²(30°) |
|---|---|---|---|
| This repository (v10.6) | 0.04 | 93.65 | 2368.93 |
| DeepSeek (rev 43) | 0.0403 | 93.74 | 2369.45 |
| Qwen (rev 44) | 0.04028 | 93.7396 | 2369.4494 |
| GLM (rev 46) | 0.040282 | 93.739590 | 2369.449491 |

The 0.1%-level offsets are convention-level (minimiser, quadrature details);
one cross-check (rev 43) localised a genuine bug in our sound-speed formula
(fixed as FIX-24). Verification artefacts live in `verification/`.

## Data

- DESI DR2 BAO: embedded (13 points + covariance, from
  [CobayaSampler/bao_data](https://github.com/CobayaSampler/bao_data), DESI
  DR2, arXiv:2503.14738).
- Pantheon+: downloaded by `data/download_data.sh` from the public
  [DataRelease](https://github.com/PantheonPlusSH0ES/DataRelease)
  (Brout et al. 2022; Scolnic et al. 2022). Not committed to the repo.

## Citation

See `CITATION.cff`. If you use this code, please cite the paper and this
repository (a Zenodo DOI will be attached to releases).

## License

MIT (code). The Pantheon+ data remain under their original release terms.
