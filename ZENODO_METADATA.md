# Zenodo metadata for the v28.0 public working draft

Use a **new version** of the existing Zenodo record. Do not overwrite or delete
the published v10.6 files.

- **Title:** Per-Probe Decomposition of Bimetric Cosmology: BAO Blindness and a Baseline Scalar-Sector Diagnostic
- **Version:** v28.0-public-draft.1
- **Publication date:** 2026-11-01
- **Author:** Aleksandr Chernikov
- **Resource type:** Software / publication package (retain the existing record's resource type if Zenodo does not permit changing it)
- **License:** MIT for code; external data remain under their original terms
- **Record-family DOI:** 10.5281/zenodo.23292159
- **Keywords:** bimetric gravity; BAO; DESI DR2; Pantheon+; compressed CMB; profile likelihood; reproducibility; working draft

## Description

This release contains manuscript v28.0, a public working draft of the
bimetric-cosmology per-probe study, together with versioned analysis archives
and historical verification materials. The methodological audit is still in
progress; this is not a final model-selection paper and does not claim an
accepted production P3/MCMC posterior.

The draft distinguishes profiling from marginalisation, limits per-probe
p-values and arithmetic χ² sums to descriptive diagnostics, and treats the
compressed-CMB result as a conditional GR-calibrated geometric diagnostic
rather than a full Boltzmann/ISW test. Baseline scalar-sector and Ω_m
statements are restricted to the tested approximations and parameter domains.
The earlier unsupported quantitative late-ISW/S₈ bound Δχ²_CMB ≲ 10⁻³ is
withdrawn.

Executable archive v10.6 remains the legacy companion to manuscripts
v21.5/v21.6. Archive v13.3 is supplied as a documentary source archive. This
release does not claim that either archive reproduces every number in v28.0
from scratch. See REPRODUCIBILITY.md for exact version boundaries.

## Recommended related identifiers

- GitHub repository: https://github.com/chernicovalex-prog/bao-blindness-bimetric
- Earlier Zenodo version / record family: https://doi.org/10.5281/zenodo.23292159

## Upload check

Include the new manuscript (Markdown and PDF), README, CHANGELOG,
REPRODUCIBILITY, release notes, citation metadata, documentary archive v13.3,
available figures, and the preserved legacy repository files. Do not include
`.git`, local credentials, internal audit correspondence, or temporary files.
