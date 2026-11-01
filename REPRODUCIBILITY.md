# Reproducibility and version boundaries

## Purpose

This document distinguishes executable reproduction of a historical release
from documentary preservation of later work. A manuscript version and a code
archive version must always be cited together.

## Version matrix

| Manuscript | Code / archive | What is claimed |
|---|---|---|
| v21.5 / v21.6 | executable v10.6 | Legacy pipeline and historical regression numbers; this is the original executable release pair |
| v28.0 public working draft | v10.6 retained; v13.3 documentary archive included | Current scientific account with audit-qualified conclusions; **not** full one-command numerical reproduction |

The files of earlier releases are retained rather than overwritten. The DOI
`10.5281/zenodo.23292159` identifies the Zenodo record family; individual
versions may also receive version-specific DOIs from Zenodo.

## Reproducing the legacy v10.6 pipeline

Environment: Python with the dependencies in `requirements.txt`.

```bash
pip install -r requirements.txt
make test
make bao
make cmb
bash data/download_data.sh && make sne-full
```

These commands concern the v10.6/v21.x release. The file
`canonical_numbers.json` contains the corresponding historical regression
values and tolerances. Pantheon+ data are external and must be obtained from
the official release through `data/download_data.sh`.

The legacy pipeline includes the massive-neutrino NR/UR transition at
`1+z = 131`. That convention belongs to the legacy release. It is not the
exact-Fermi–Dirac P3 target now under methodological audit.

## What v13.3 provides

`Bimetric_Code_Archive_v13.3.txt` preserves the later source state for
inspection and provenance. It is supplied as a documentary text archive.
This release does not claim that:

- it is packaged as a clean, standalone executable environment;
- every later diagnostic has an independently validated production runner;
- it reproduces every edited number or statement in v28.0;
- it supplies an accepted production P3/MCMC posterior;
- it provides a full Boltzmann, recombination, or late-ISW calculation.

Accordingly, the README no longer says that this repository reproduces every
number of the current manuscript from scratch.

## Evidence boundaries in manuscript v28.0

1. **Profiling is not marginalisation.** The per-probe `f` and `M` treatment
   produces profile-likelihood diagnostics. It does not define an integrated
   marginal posterior, evidence, or Bayes factor.
2. **Per-probe sums are descriptive.** Adding reported χ² contributions is
   arithmetic bookkeeping unless a common, audited joint likelihood and
   nuisance treatment are supplied.
3. **Compressed CMB is conditional.** The calculation uses a GR-calibrated
   geometric diagnostic. Its physical transfer error for the bimetric model
   has not been quantified by an end-to-end Boltzmann/ISW computation.
4. **No accepted production P3 posterior is claimed.** The current public
   draft must not be read as a completed MCMC/model-selection analysis.
5. **Baseline exclusions are bounded.** Statements about the BAO baseline,
   θ-grid, and Ω_m rigidity refer to the tested conventions, grids, and local
   analytic approximations, not to the entire bimetric parameter space.
6. **The late-ISW/S₈ numerical bound is withdrawn.** The former
   `Δχ²_CMB ≲ 10⁻³` estimate lacked an end-to-end justification and is not a
   reproducibility target.

## Figures

The package includes the available figure files `Fig3_theta_scan.png` and
`FigD1_SNe_reference_sensitivity.png`. Captions for figures whose original
image files are not present in this package are retained in the manuscript;
no replacement scientific data have been generated for this editorial
release.

## Recommended citation practice

Cite all applicable identifiers explicitly, for example:

- scientific claim: manuscript v28.0, public working draft, 1 November 2026;
- legacy execution: code archive v10.6 and its corresponding v21.x manuscript;
- later source inspection: documentary archive v13.3;
- persistent archive: the relevant Zenodo version DOI, or the record-family
  DOI `10.5281/zenodo.23292159` when a version DOI is unavailable.
