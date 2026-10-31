# BAO Blindness to Bimetric Gravity

Public working draft and versioned analysis archive for the bimetric-cosmology
per-probe study by Aleksandr Chernikov.

**Current manuscript:** [v28.0, public working draft](BAO_Blindness_v28.0_public_draft.md),
released 1 November 2026. The methodological audit is still in progress.

## Read this first

This repository contains materials from different stages of one continuing
research programme. They are **not interchangeable**:

- **v28.0 public working draft** is the current public account of the study.
  It is not a final model-selection paper.
- **Code archive v10.6** is the executable companion to manuscripts
  v21.5/v21.6. It remains available as a legacy pipeline.
- **Code archive v13.3** is included as a documentary source archive. It is
  not claimed to be a drop-in executable that reproduces every number in
  v28.0.
- Manuscripts v21.5 and v21.6 and archive v10.6 are retained for provenance;
  they are not silently rewritten.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the exact version boundaries
and [CHANGELOG.md](CHANGELOG.md) for the public release history.

## Scientific status of v28.0

The draft decomposes bimetric-cosmology diagnostics by probe: DESI DR2 BAO,
Pantheon+ SNe Ia, and a compressed-CMB geometric diagnostic. Its main
conclusions are deliberately limited:

- the reported BAO result concerns the tested baseline and per-probe diagnostic;
- the compressed-CMB calculation is a conditional geometric diagnostic, not a
  full Boltzmann/ISW calculation;
- `f` and `M` are profiled, not marginalised;
- per-probe profile likelihoods are not an integrated marginal posterior,
  evidence, or Bayes factor;
- descriptive arithmetic sums of per-probe χ² contributions are not a
  certified joint likelihood or a model-selection p-value;
- no accepted production P3/MCMC posterior is claimed;
- a matched model comparison is still being audited.

The earlier quantitative statement that late-ISW/S₈ effects must satisfy
`Δχ²_CMB ≲ 10⁻³` is not supported by an end-to-end calculation and is
explicitly withdrawn in v28.0.

## Repository layout

| Path | Content and status |
|---|---|
| `BAO_Blindness_v28.0_public_draft.md` | Current public working draft |
| `BAO_Blindness_v28.0_public_draft.pdf` | PDF rendering of the same draft |
| `BAO_Blindness_v21.5.md`, `BAO_Blindness_v21.6.md` | Historical manuscripts |
| `Bimetric_Code_Archive_v10.6.txt` | Legacy executable source archive |
| `Bimetric_Code_Archive_v13.3.txt` | Documentary source archive |
| `bimetric/`, `probes/`, `tests/` | Legacy v10.6 pipeline |
| `canonical_numbers.json` | Canonical numbers for the legacy v10.6/v21.x pair |
| `verification/` | Historical independent-verification artefacts |
| `figures/` | Available manuscript figures |
| `REPRODUCIBILITY.md` | Version-specific reproduction instructions |

## Legacy quick start — v10.6 only

The following commands apply to the legacy v10.6 pipeline, **not to a full
reproduction of manuscript v28.0**:

```bash
pip install -r requirements.txt
make test        # legacy self-checks and v10.6 regressions
make bao         # legacy DESI DR2 BAO diagnostic
make cmb         # legacy compressed-CMB theta scan
bash data/download_data.sh && make sne-full
```

The legacy background implementation uses the massive-neutrino transition at
`1+z = 131`. This convention must not be confused with the exact-Fermi–Dirac
P3 target discussed in the current draft and its ongoing audit.

Pantheon+ data are downloaded from the official public release and are not
committed here. Their original release terms continue to apply.

## Historical canonical results — archive v10.6

| Probe | Quantity | Legacy value |
|---|---|---|
| Baseline | y(0), Ω_m, E²(0)−1 | 0.999512, 0.30439, 5·10⁻¹⁵ |
| BAO | Δχ², full covariance | 0.178 |
| SNe | Δχ², diagonal / full | 15.58 / 6.50 |
| CMB | Δχ² at θ = 6.1° / 30° | 0.04 / 2368.93 |

These values belong to the legacy v10.6/v21.x release. They are retained for
provenance and regression testing, not presented as a completed v28.0
production analysis.

## Verification history

The legacy CMB θ-scan was checked by independent codebases written from a
blind specification and executed by the project administrator. The historical
artefacts and logs are in `verification/`. Independent implementation checks
support particular numerical elements; they do not by themselves validate the
physical transfer of a GR-calibrated compressed likelihood or establish a
full model-selection result.

## Citation and archive

See [CITATION.cff](CITATION.cff). The versioned Zenodo record is
[10.5281/zenodo.23292160](https://doi.org/10.5281/zenodo.23292160).
Please cite the manuscript version and code/archive version actually used.

## License

MIT for the code, unless a file states otherwise. External data remain under
their original release terms.
