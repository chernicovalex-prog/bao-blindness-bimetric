# Dependency map for the corrected candidate

This map records which artifacts must move together. It supersedes the pre-integration numeric dependency map; the older article/code and their numeric values are not mixed into this candidate.

## Main dependency chain

```text
DESI DR2 BAO + Pantheon+ + Bansal–Huterer CMB + H&M benchmark
  -> core/probe backgrounds
  -> CMB/BAO/SNe raw fits and scans
  -> summary.json and per-probe tables
  -> anchor_scan.json and scalar_at_minima.json
  -> exact-stream MC and cmb_linbias_check.json
  -> canonical assertions (checks.py)
  -> figures and article numbers
  -> README, reproduction commands and article build
```

## Inputs and controls

| Input/convention | Downstream rule |
|---|---|
| DESI DR2 BAO mean/covariance | BAO baseline/scans, anchor BAO column and shared comparator. Full published null is `synthetic_null_refit(n=2000, seed=20261014)`. |
| Pantheon+ curated fields/covariance | SNe baseline/scans, anchor SNe column and every MC noise realisation. The curated header is not an official byte-identical claim. |
| Bansal–Huterer vector/covariance | CMB baseline, theta/mFP/(alpha,beta) grids, lambda sweeps, anchor CMB column and CMB linearisation. |
| `omega_ng(z, Om, h)` | Must be identical in quartic, E2, z=0 closure, direct/spline, split and live CMB comparators. Cache identity includes `h`. |
| Hu–Sugiyama 39.5 | Recombination redshift for both bimetric and LCDM CMB paths. |
| BAO photons-only vs SNe/scalar fixed radiation | Keep separate regression controls; do not propagate the CMB h-dependent density into BAO/SNe/scalar fixed conventions. |
| Reduced $\varepsilon=0$ certificate vs full-model scan | Exact reduced algebraic proof and full-model numerical evidence must not be conflated. |

## Raw artifact dependencies

| Upstream change | Artifacts requiring coordinated regeneration |
|---|---|
| SCI-01 density + SCI-10 coefficient | `summary.json`, `cmb_grids.json`, both lambda sweeps, split results, anchor tree, lambda=1 profile, shared comparator and every CMB-dependent article value. |
| Corrected CMB + unified statistic | Fine/offset calibration, per-realisation arrays, coverage, argmin, dchi_true and exact-vs-linear bias. Old coverage is not inherited. |
| Changed anchor/profile minima | Continuous rows, per-probe decomposition, `scalar_at_minima.json`, anchor tables/figures and all repeated article values. |
| Changed BAO/SNe/scalar fixed inputs | Regression controls and dependent tables; CMB-only changes must not alter them. |
| Pantheon header only | Provenance text only; all 1701 data rows and covariance hashes remain controls. |

## MC exact-stream dependencies

- F0–F3 are contiguous segments of one random stream with `seeds=[20261118]`.
- Each segment checkpoint stores `done`, per-realisation arrays and RNG-state.
- Central merge depends on all four complete segments and identical code/data signatures.
- `profile_calibration_offset.json` is independent of the fine merge but uses the same unified statistic and covariance convention.
- `cmb_linbias_check.json` is independent exact-vs-linear evidence and does not convert the linearised path into exact MC.

## Deferred branches

- `joint_tie.py` and `joint_tie_l1.py`: `DEFERRED_NO_LIVE_ARTICLE`. The live shared comparator is `anchor_scan.json:lcdm_shared`, not a tied-rd output.
- `joint_tie_hm.py` and `joint_tie_hm_ext.py`: deferred H&M reconstruction, no runs and no imported numbers.
- Improved compressed-CMB method, Boltzmann code, MCMC and optional robustness surveys are outside the current integration.

## Text and documentation dependencies

`G-06/G-07` artifacts feed `checks.py` expectations, figure assertions and article values. G-09 synchronizes assertions, G-10 applies approved SCI text and bibliography, G-11 verifies code/figures, G-12 checks editorial consistency, G-13 verifies the 35-page build and G-14 freezes these documents. A later scientific REV must use the complete dependency tree rather than selected numbers.
