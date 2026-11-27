# The H&M Anchor, Decomposed — current research preview

This directory contains the current laboratory candidate of the per-probe bimetric analysis by Aleksandr Chernikov.

**Status:** research preview, not peer reviewed. No final scientific acceptance or full CMB validation is claimed.

## Entry points

- `PRERELEASE_README.md` — package contents and verification summary.
- `src/` — current Python implementation.
- `results/` — declared current results used by the assertions and figures.
- `data/` — small public DESI inputs and `download_data.sh`, which retrieves the official Pantheon+ table and pins the reviewed covariance by SHA-256.
- `01_Statya/` — article source; figures are regenerated from declared results by CI.
- `package_release.py` — builds the PDF, arXiv source, source/data package, checksums and full research-preview ZIP.

## Reproduction

Start with `README.md` and `requirements.txt`. The documented route is staged and checkpointed. Every stage has its own RC, logs and evidence; no uninterrupted monolithic production run is claimed.

The CI workflow verifies the CMB gates, independent background, current assertions and figures, then builds the article and packages. A GitHub pre-release is created only by an explicit manual workflow dispatch with `release=true`.

## Scientific boundary

The work separates background fit, viability and model preference. Compressed CMB remains a conditional GR-calibrated diagnostic, not a full perturbation-level bimetric CMB test. The H&M joint MAP values and the descriptive per-probe sum are different estimators and are not interchangeable.

Earlier versions are preserved under `archive/`. In particular, articles v21.5/v21.6 belong to code v10.6, while article v28.0/v28.0.1 belongs to documentary code v13.3. These version pairs must not be mixed.

## License and data

MIT for code. DESI DR2 BAO, Pantheon+ and Bansal–Huterer compressed-CMB inputs are public; see `DATA_PROVENANCE.md`. No secrets, credentials, internal journals or blind-isolation material are included.
