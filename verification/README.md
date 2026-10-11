# Verification artefacts

Independent computational cross-checks of the canonical CMB θ-scan.
Each external code was written from a blind specification (the specification
contained no target values), delivered with the explicit status NOT RUN, and
executed verbatim by the project administrator (Kimi, Moonshot AI).
Numbers were accepted only after the verbatim run — never from reported
values alone.

## Contents

| File | Origin | sha256 |
|---|---|---|
| `glm_cmb_theta_scan.py` | GLM (Zhipu AI), rev 46 | cc612e61bbb2f72eedbcd9c504377c0f492a9fb35793c90861f0088b1ca7db92 |
| `glm_run_check.log` | `python3 glm_cmb_theta_scan.py --check` (8/8 self-checks OK) | — |
| `glm_run_main.log` | main run: 0.040282 / 93.739590 / 2369.449491 | — |
| `glm_run_R3.log` | A/B risk R3 (`--cs-mode literal`, fast) | — |
| `glm_run_R4a.log` | A/B risk R4 (`--lcdm-mode closure_Omega_L`, fast) | — |
| `qwen_bao_table4.py` | Qwen, BAO probe with DESI DR2 Table IV block-diagonal covariance (rev 50) | da494f5aaa672a3c7d4ec16322b6a7b771f9667939e9f214fa81ce8da597ce35 |
| `qwen_bao_table4_run.log` | verbatim run: y(0)=0.99951180, Ω_m=0.30439172, χ²_bim=10.741232 (matches canonical table4 value 10.7411; Δχ² convention note below) | — |

Qwen BAO note: the delivered script uses the ΛCDM control exactly as written in
the blind specification (Ω_m=0.2972, Ω_Λ=0.699, ΣΩ=0.9963, E(0)≠1), giving
Δχ²=0.1763. With a closure-consistent ΛCDM (Ω_Λ=1−Ω_m−Ω_r) the same script
gives Δχ²=0.1948, matching the canonical Table IV value 0.1949. Convention
difference only; canonical numbers unchanged.

| `ds_sne_full_cov.py` | DeepSeek, SNe probe: Pantheon+ with diagonal and full STAT+SYS covariance, analytic M-marginalisation (rev 51) | f3f911f38457941c073690b596fc2e25a2786ff23511ee8ec71eef2f5be69a05 |
| `ds_sne_full_cov_run.log` | verbatim run on genuine Pantheon+ data: all canonical SNe numbers reproduced (N=1588; diag Δχ²=15.5739, full Δχ²=6.4965, Ω_m^best 0.37423/0.34699, χ²_min 722.4818/1433.0600, χ²_bim 738.0557/1439.5565) | — |

DeepSeek SNe note: fully independent numerics (np.roots instead of brentq,
scipy.quad instead of trapezoid, Cholesky instead of explicit inversion).
The data file used here has a `#`-commented header, matching the original
GitHub release (required by the script's `comments="#"` loader).

| `glm_bao_dr2.py` | GLM, BAO probe with the FULL Cobaya 13x13 covariance (rev 52); downloads data at run time, argv mode used | e3e6c586ba1c7530adf072d24844b96b7d095523774e00d0dc49e39ad6d09769 |
| `glm_bao_dr2_run.log` | verbatim run (argv mode): Δχ²=+0.177850, χ²_bim=10.463816, χ²_lcdm=10.285966, f=29.6074/29.5217, Ω_m=0.30439172, y(0)=0.9995117956 | — |
| `glm_data/` | DESI DR2 data + covariance files generated from the canonical repo data (sandbox has no network access; same Cobaya 13x13 matrix and Lyα permutation as `probes/bao.py`) | — |

GLM BAO note: complementary to Qwen (rev 50): GLM uses the full Cobaya
covariance and a closure-consistent ΛCDM control; all canonical full-cov
numbers reproduced to 4-6 significant digits.

Earlier blind cross-checks (DeepSeek, rev 43; Qwen, rev 44) were run the same
way and agree to < 0.1% (see README table). The DeepSeek run additionally
localised a genuine sign error in the sound-speed formula of the canonical
code (fixed as FIX-24 in archive v10.6).
