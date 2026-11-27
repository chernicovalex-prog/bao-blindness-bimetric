# Bimetric cosmology: per-probe analysis and research previews

Public repository for Aleksandr Chernikov's bimetric-cosmology work.

## Current research preview

The current work is **The H&M Anchor, Decomposed: A Per-Probe Analysis of the Bimetric Cosmological Benchmark**.

See [`current/README.md`](current/README.md) and [`current/PRERELEASE_README.md`](current/PRERELEASE_README.md).

The manuscript decomposes the Högås–Mörtsell bimetric benchmark across DESI DR2 BAO, Pantheon+ supernovae and a conditional compressed-CMB diagnostic, and evaluates the scalar gradient-instability boundary. A good background fit is not treated as proof of viability or model preference.

This is a research preview, not a peer-reviewed paper. The GR-calibrated compressed-CMB transfer remains conditional; full perturbation-level CMB validity and a full H&M-convention joint reconstruction are not established here.

## Repository layout

- `current/` — current sources, declared results, article inputs, data downloader and packaging workflow.
- `archive/` — historical articles, code, figures and logs, retained without overwriting.
- `bimetric/`, `probes/`, `tests/`, `verification/` — earlier executable and verification lines, retained for their documented legacy scope.

## Version boundaries

- Articles **v21.5/v21.6** belong to executable code **v10.6**.
- Article **v28.0/v28.0.1** belongs to documentary code **v13.3**.
- These pairs are not interchangeable and must not be mixed.

The earlier executable v10.6 pipeline remains a legacy line. Archive v13.3 is documentary and is not claimed to reproduce every number of the current candidate without a separate check.

## Citation and archive

The existing Zenodo concept DOI is [10.5281/zenodo.23292159](https://doi.org/10.5281/zenodo.23292159). When a new research-preview version is published, use the version-specific DOI displayed by Zenodo. Legacy executable v10.6 has the version-specific DOI `10.5281/zenodo.23292160`.

## License

MIT for code. External public data retain their original terms. No credentials, internal journals or blind-isolation material belong to the current release package.
