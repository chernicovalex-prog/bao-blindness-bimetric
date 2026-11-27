# Data provenance

- DESI DR2 BAO: the 13 distance measurements of DESI DR2 Results II Table IV and the full 13x13 covariance as distributed by CobayaSampler `bao_data` (`desi_bao_dr2`). Public.
- Pantheon+: the official distance table and STAT+SYS covariance from the PantheonPlusSH0ES DataRelease are downloaded by `data/download_data.sh`. Only the table's column-name line is normalized to a comment; all 1701 data rows remain unchanged. The pipeline reads `zHD`, `zCMB`, `zHEL`, `m_b_corr` and `m_b_corr_err_DIAG`. The covariance SHA-256 is pinned to the reviewed public file; `src/checks.py` verifies all quoted SNe results and event counts.
- Compressed CMB: Bansal–Huterer vector and covariance embedded in `src/probes/cmb.py` (published values). CMBComp is discussed but not substituted.
- Hassan–Rosen benchmark and literature references are cited in `01_Statya/refs.bib`.

No secrets, access tokens, cookies or private credentials are included.
