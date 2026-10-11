#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TASK 1 (BAO, DESI DR2): Delta chi2 = chi2(bimetric) - chi2(LCDM)
================================================================
STATUS: NOT RUN.
  The author has no execution environment; this file is handed over for
  verbatim execution.  It contains NO transcribed DESI values and NO
  pre-filled results: every measured number enters at run time from the
  official data source, so nothing can be imitated.

MODEL (bimetric background; Hogas & Mortsell 2025, arXiv:2507.03743)
  Flat background; y = a_f/a_g = physical root in (0,1] of quartic (A13),
  with m2 = (m_FP/H0)^2, c2 = cos^2(theta), s2 = sin^2(theta):
    a0 = -(1/3) c2 m2 (1 + 2 alpha + beta)
    a1 = Om_tot + Om_L + m2 [ c2 (alpha + beta) - s2 (1 + alpha + beta/3) ]
    a2 = m2 [ -c2 beta + s2 (1 + 2 alpha + beta) ]
    a3 = -[ Om_L + (1/3) m2 ( c2 (-1 + alpha - beta) + 3 s2 (alpha + beta) ) ]
    a4 = (1/3) s2 m2 beta
    P(y) = a0 + a1 y + a2 y^2 + a3 y^3 + a4 y^4 = 0
    Om_tot(z) = Om_m (1+z)^3 + Om_r (1+z)^4
    Om_DE(z)  = Om_L - s2 m2 (1-y) [ 1 + alpha (1-y) + (beta/3) (1-y)^2 ]
    E^2(z)    = Om_tot(z) + Om_DE(z)
  Fixed by the task: theta = 6.1 deg, m_FP/H0 = 25, alpha = 0, beta = 3,
  Om_L = 0.699, Om_r = 5.38e-5; Om_m from the closure E^2(0) = 1.
  LCDM control: flat, Om_m = 0.2972 (fixed), Om_L = 1 - Om_m - Om_r.

DATA (choice stated, as required)
  DESI DR2 BAO (arXiv:2503.14738, Table IV): 13 points = BGS D_V +
  (D_M, D_H) pairs for LRG1, LRG2, LRG3+ELG1, ELG2, QSO, Lya.
  COVARIANCE CHOICE: full 13x13 matrix from the repository
  CobayaSampler/bao_data (desi_bao_dr2), auto-discovered and downloaded
  at run time.  Also accepted: matching local files in the working
  directory, or explicit sources:
      python this_script.py <data_txt_path_or_url> <cov_path_or_url>
  The alternative (block-diagonal covariance hand-built from Table IV
  sigmas/correlations) is deliberately NOT used, to avoid transcribing
  measured numbers by hand.

FREE PARAMETER f = c/(H0 r_d), MARGINALIZED BY MINIMIZATION
  Predictions are exactly linear in f:
    D_M/r_d = f I(z),  I(z) = Int_0^z dz'/E(z')
    D_H/r_d = f / E(z)
    D_V/r_d = f ( z I(z)^2 / E(z) )^(1/3)
  hence the profile minimum is analytic:
    f_hat    = (v' C^-1 d) / (v' C^-1 v)
    chi2_min = d' C^-1 d - (v' C^-1 d)^2 / (v' C^-1 v)
    sigma_f  = (v' C^-1 v)^(-1/2)
  A scipy minimize_scalar cross-check is printed.

PRINTED OUTPUT
  Om_m (closure), y(0), chi2_bimetric, chi2_LCDM,
  Delta chi2 = chi2_bimetric - chi2_LCDM, best f (both models),
  chi2/dof (dof = 13 - 1), full parsed-data table, consistency checks
  (E^2(0)-1; analytic vs numerical f-minimization; txt errors vs
  covariance diagonal when the data file provides errors).

DEPENDENCIES: python3, numpy, scipy.  Nothing else.
"""

import json
import os
import re
import sys
import urllib.request

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

# ----------------------------------------------------------------------
# 1. Fixed parameters (task specification)
# ----------------------------------------------------------------------
THETA_DEG = 6.1
M_FP_OVER_H0 = 25.0
ALPHA = 0.0
BETA = 3.0
OM_LAM = 0.699
OM_R = 5.38e-5
OM_M_LCDM = 0.2972

M2 = M_FP_OVER_H0 ** 2
_TH = np.radians(THETA_DEG)
C2 = float(np.cos(_TH) ** 2)
S2 = float(np.sin(_TH) ** 2)

# ----------------------------------------------------------------------
# 2. Bimetric background: quartic (A13) -> y(z), Om_DE(z), E^2(z)
# ----------------------------------------------------------------------

def quartic_coeffs(Om_tot):
    """Coefficients a0..a4 of P(y) = a0 + a1 y + a2 y^2 + a3 y^3 + a4 y^4."""
    a0 = -(1.0 / 3.0) * C2 * M2 * (1.0 + 2.0 * ALPHA + BETA)
    a1 = (Om_tot + OM_LAM
          + M2 * (C2 * (ALPHA + BETA) - S2 * (1.0 + ALPHA + BETA / 3.0)))
    a2 = M2 * (-C2 * BETA + S2 * (1.0 + 2.0 * ALPHA + BETA))
    a3 = -(OM_LAM + (1.0 / 3.0) * M2
           * (C2 * (-1.0 + ALPHA - BETA) + 3.0 * S2 * (ALPHA + BETA)))
    a4 = (1.0 / 3.0) * S2 * M2 * BETA
    return a0, a1, a2, a3, a4


def y_physical(Om_tot):
    """Physical root y in (0,1].

    Exact identities for this coefficient set: P(0) = a0 < 0 and
    P(1) = Om_tot > 0 (the m2-terms cancel identically at y = 1), and
    P'(y) > 0 on [0,1] for all Om_tot >= 0 with the given parameters,
    so the root in (0,1) exists and is unique; brentq needs no guess.
    """
    a0, a1, a2, a3, a4 = quartic_coeffs(Om_tot)

    def P(y):
        return a0 + y * (a1 + y * (a2 + y * (a3 + y * a4)))

    if not (P(0.0) < 0.0 < P(1.0)):
        raise RuntimeError("quartic bracket failed for Om_tot=%r" % (Om_tot,))
    return brentq(P, 0.0, 1.0, xtol=1e-15, rtol=1e-15, maxiter=300)


def Om_DE_of_y(y):
    """Effective dark-energy density parameter for given y."""
    d = 1.0 - y
    return OM_LAM - S2 * M2 * d * (1.0 + ALPHA * d + (BETA / 3.0) * d * d)


def make_E2_bim(Om_m):
    """E^2(z) of the bimetric background for given Om_m."""
    def E2(z):
        zp1 = 1.0 + z
        Om_tot = Om_m * zp1 ** 3 + OM_R * zp1 ** 4
        return Om_tot + Om_DE_of_y(y_physical(Om_tot))
    return E2


def E2_lcdm(z):
    """E^2(z) of the flat LCDM control (Om_m = 0.2972 fixed)."""
    zp1 = 1.0 + z
    return (OM_M_LCDM * zp1 ** 3 + OM_R * zp1 ** 4
            + (1.0 - OM_M_LCDM - OM_R))


def solve_Om_m_closure():
    """Bimetric Om_m from the closure E^2(0) = 1 (all else fixed)."""
    def resid(Om_m):
        return (Om_m + OM_R
                + Om_DE_of_y(y_physical(Om_m + OM_R)) - 1.0)

    lo, hi = 1e-8, 0.9
    if resid(lo) * resid(hi) > 0.0:
        raise RuntimeError("no bracket for the closure on [%g, %g]" % (lo, hi))
    return brentq(resid, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=300)


# ----------------------------------------------------------------------
# 3. Theory vectors and profile likelihood in f = c/(H0 r_d)
# ----------------------------------------------------------------------

def theory_vector(E2, zs, obs):
    """v such that prediction_i = f * v_i, with f = c/(H0 r_d)."""
    v = np.empty(len(zs), dtype=float)
    for i in range(len(zs)):
        z = float(zs[i])
        invE = 1.0 / np.sqrt(E2(z))
        if obs[i] == "DH":
            v[i] = invE
        elif obs[i] == "DM":
            v[i] = quad(lambda t: 1.0 / np.sqrt(E2(t)), 0.0, z,
                        epsabs=1e-13, epsrel=1e-11, limit=200)[0]
        elif obs[i] == "DV":
            I = quad(lambda t: 1.0 / np.sqrt(E2(t)), 0.0, z,
                     epsabs=1e-13, epsrel=1e-11, limit=200)[0]
            v[i] = (I * I * z * invE) ** (1.0 / 3.0)
        else:
            raise ValueError("unknown observable %r" % (obs[i],))
    return v


def profile_f(d, C, v):
    """Analytic minimum of chi2(f) = (d - f v)' C^-1 (d - f v)."""
    Wd = np.linalg.solve(C, d)
    Wv = np.linalg.solve(C, v)
    vWv = float(np.dot(v, Wv))
    vWd = float(np.dot(v, Wd))
    dWd = float(np.dot(d, Wd))
    f_hat = vWd / vWv
    return dWd - vWd * f_hat, f_hat, vWv ** -0.5


def chi2_at_f(f, d, C, v):
    r = d - f * v
    return float(np.dot(r, np.linalg.solve(C, r)))


# ----------------------------------------------------------------------
# 4. Data acquisition (no transcribed numbers anywhere)
# ----------------------------------------------------------------------
REPO = "CobayaSampler/bao_data"
API = "https://api.github.com/repos/%s/contents" % REPO
RAW_BASES = ["https://raw.githubusercontent.com/%s/%s/" % (REPO, br)
             for br in ("main", "master")]
RAW_STEMS = [
    "desi_bao_dr2/desi_bao_dr2",
    "desi_dr2/desi_dr2_bao",
    "desi_dr2/desi_dr2",
    "desi_bao_dr2/desi_bao_dr2_bao",
    "desi_dr2/desi_bao_dr2",
    "desi_bao_dr2_bao",
    "desi_dr2_bao",
    "desi_bao_dr2",
]
LOCAL_TXT = ["desi_bao_dr2.txt", "desi_dr2_bao.txt", "desi_dr2.txt",
             "desi_bao_dr2_bao.txt"]
LOCAL_COV = [s[:-4] + ".cov" for s in LOCAL_TXT]

_NUM = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


def _is_num(tok):
    return _NUM.match(tok) is not None


def _classify(tok):
    """Map an observable token to 'DM' | 'DH' | 'DV' (or None)."""
    u = tok.upper().replace("_", "").replace("-", "").replace("/", "")
    for key in ("DM", "DH", "DV"):
        if key in u and ("RD" in u or "OVER" in u or u == key):
            return key
    return None


def parse_data_text(txt):
    """Parse data lines 'z_eff observable value [error]' (order-tolerant).

    Returns a list of (z_eff, 'DM'|'DH'|'DV', value, error-or-None).
    """
    rows = []
    for line in txt.splitlines():
        s = line.strip()
        if not s or s[0] in "#%!;":
            continue
        toks = s.split()
        i_obs = None
        for i, t in enumerate(toks):
            if not _is_num(t):
                k = _classify(t)
                if k is not None:
                    i_obs = i
                    break
        if i_obs is None:
            continue
        before = [t for t in toks[:i_obs] if _is_num(t)]
        after = [t for t in toks[i_obs + 1:] if _is_num(t)]
        z = val = err = None
        if before and after:
            z, val = float(before[-1]), float(after[0])
            err = float(after[1]) if len(after) > 1 else None
            if not (0.0 < z < 5.0) and len(before) >= 2 \
                    and 0.0 < float(before[-2]) < 5.0:
                # format 'z value name [err]': z and value precede the name
                z, val = float(before[-2]), float(before[-1])
                err = float(after[0])
        elif len(after) >= 2:
            z, val = float(after[0]), float(after[1])
            err = float(after[2]) if len(after) > 2 else None
        elif len(before) >= 2:
            z, val = float(before[-2]), float(before[-1])
        if z is None or val is None:
            continue
        rows.append((z, _classify(toks[i_obs]), val, err))
    return rows


def parse_cov_text(txt, n):
    """Parse a covariance file: full n x n matrix, or a triangle, with an
    optional leading dimension token.  Positive-definiteness is checked."""
    nums = []
    for line in txt.splitlines():
        s = line.strip()
        if not s or s[0] in "#%!;":
            continue
        for tok in s.replace(",", " ").split():
            if _is_num(tok):
                nums.append(float(tok))
    nums = np.asarray(nums, dtype=float)
    full, tri = n * n, n * (n + 1) // 2
    if nums.size not in (full, tri):
        if nums.size == full + 1 and abs(nums[0] - n) < 1e-12:
            nums = nums[1:]
        elif nums.size == tri + 1 and abs(nums[0] - n) < 1e-12:
            nums = nums[1:]
    if nums.size == full:
        M = nums.reshape(n, n)
        C = 0.5 * (M + M.T)
    elif nums.size == tri:
        C = None
        for mode in ("lower", "upper"):
            M = np.zeros((n, n))
            k = 0
            if mode == "lower":
                for i in range(n):
                    for j in range(i + 1):
                        M[i, j] = nums[k]
                        k += 1
            else:
                for j in range(n):
                    for i in range(j + 1):
                        M[i, j] = nums[k]
                        k += 1
            M = 0.5 * (M + M.T)
            try:
                np.linalg.cholesky(M)
                C = M
                break
            except np.linalg.LinAlgError:
                pass
        if C is None:
            raise ValueError("triangular cov parsed but not positive definite")
    else:
        raise ValueError("cov file has %d numbers; expected %d (full %dx%d) "
                         "or %d (triangle)" % (nums.size, full, n, n, tri))
    try:
        np.linalg.cholesky(C)
    except np.linalg.LinAlgError:
        raise ValueError("parsed covariance is not positive definite")
    return C


def fetch_text(src):
    """Read a local file or download a URL; return text."""
    if os.path.exists(src):
        with open(src, "r", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    req = urllib.request.Request(src, headers={"User-Agent": "bao-dchi2"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")


def candidate_urls():
    """Candidate (txt, cov) URLs: GitHub contents API discovery first
    (adapts to the actual naming in CobayaSampler/bao_data), then fixed
    raw-URL guesses on both common branch names."""
    txts, covs = [], []
    try:
        top = json.loads(fetch_text(API))
        for e in (top if isinstance(top, list) else []):
            if not isinstance(e, dict):
                continue
            name = str(e.get("name", ""))
            path = str(e.get("path", ""))
            if e.get("type") == "dir" and re.search(r"desi", name, re.I):
                dir_is_dr2 = bool(re.search(r"dr2", name, re.I))
                try:
                    sub = json.loads(fetch_text(API + "/" + path))
                except Exception:
                    continue
                for f in (sub if isinstance(sub, list) else []):
                    if not isinstance(f, dict) or f.get("type") != "file":
                        continue
                    fname = str(f.get("name", ""))
                    if ((dir_is_dr2 or re.search(r"dr2", fname, re.I))
                            and fname.lower().endswith((".txt", ".cov"))):
                        url = f.get("download_url") \
                            or (RAW_BASES[0] + path + "/" + fname)
                        (covs if fname.lower().endswith(".cov")
                         else txts).append(url)
            elif (e.get("type") == "file" and re.search(r"desi", name, re.I)
                  and re.search(r"dr2", name, re.I)
                  and name.lower().endswith((".txt", ".cov"))):
                url = e.get("download_url") or (RAW_BASES[0] + name)
                (covs if name.lower().endswith(".cov") else txts).append(url)
    except Exception as exc:
        print("[data] GitHub API discovery failed (%r); using fixed URL "
              "guesses" % (exc,))
    for base in RAW_BASES:
        for stem in RAW_STEMS:
            txts.append(base + stem + ".txt")
            covs.append(base + stem + ".cov")
    return txts, covs


def rows_to_table(rows):
    return "\n".join(
        "  z=%9.4f  %-2s  value=%10.5f%s"
        % (r[0], r[1], r[2], "" if r[3] is None else "  err=%.5f" % r[3])
        for r in rows)


def validate(rows, C):
    n = len(rows)
    if n != 13:
        raise ValueError("expected 13 DESI DR2 BAO points, got %d:\n%s"
                         % (n, rows_to_table(rows)))
    obs = [r[1] for r in rows]
    if (obs.count("DV"), obs.count("DM"), obs.count("DH")) != (1, 6, 6):
        raise ValueError("expected 1 DV + 6 DM + 6 DH, got DV=%d DM=%d DH=%d:"
                         "\n%s"
                         % (obs.count("DV"), obs.count("DM"), obs.count("DH"),
                            rows_to_table(rows)))
    zs = [r[0] for r in rows]
    if any((zz <= 0.0) or (zz >= 5.0) for zz in zs):
        raise ValueError("unphysical z_eff values:\n%s" % rows_to_table(rows))
    if C.shape != (13, 13):
        raise ValueError("covariance shape %r != (13, 13)" % (C.shape,))
    for i, o in enumerate(obs):
        if o in ("DM", "DH"):
            mate = "DH" if o == "DM" else "DM"
            if not any(obs[j] == mate and abs(zs[j] - zs[i]) < 5e-3
                       for j in range(n)):
                raise ValueError("no %s partner for %s at z=%.4f"
                                 % (mate, o, zs[i]))


def load_data():
    """Load the DESI DR2 BAO data (13 points) and the full 13x13 covariance.

    Priority: (1) explicit argv sources, (2) matching local files,
    (3) download from CobayaSampler/bao_data.  Returns (rows, C, source).
    """

    def structure_ok(rows):
        if len(rows) != 13:
            return False
        o = [r[1] for r in rows]
        return (o.count("DV"), o.count("DM"), o.count("DH")) == (1, 6, 6)

    def try_pair(txt, ctxt):
        rows = parse_data_text(txt)
        if not structure_ok(rows):
            return None
        try:
            C = parse_cov_text(ctxt, 13)
        except Exception:
            return None
        return rows, C

    if len(sys.argv) == 2:
        print("[data] need two arguments (data txt + cov); "
              "ignoring the single argument given")
    if len(sys.argv) >= 3:
        pair = try_pair(fetch_text(sys.argv[1]), fetch_text(sys.argv[2]))
        if pair is None:
            raise RuntimeError("argv data/cov failed the structure check "
                               "(need 13 points: 1 DV + 6 DM + 6 DH; "
                               "13x13 positive-definite cov)")
        return pair[0], pair[1], "argv: %s | %s" % (sys.argv[1], sys.argv[2])

    for t, c in zip(LOCAL_TXT, LOCAL_COV):
        if os.path.exists(t) and os.path.exists(c):
            pair = try_pair(fetch_text(t), fetch_text(c))
            if pair is not None:
                return pair[0], pair[1], "local files: %s | %s" % (t, c)
            print("[data] local pair (%s, %s) failed the structure check; "
                  "continuing to the repository" % (t, c))

    txts, covs = candidate_urls()
    tried = []
    for ut in txts:
        try:
            txt = fetch_text(ut)
        except Exception as exc:
            tried.append("%s   [%r]" % (ut, exc))
            continue
        rows = parse_data_text(txt)
        if not structure_ok(rows):
            tried.append("%s   [parsed %d data rows]" % (ut, len(rows)))
            continue
        sibling = ut.rsplit(".", 1)[0] + ".cov"
        for uc in [sibling] + [c for c in covs if c != sibling]:
            try:
                ctxt = fetch_text(uc)
                C = parse_cov_text(ctxt, 13)
            except Exception:
                continue
            return rows, C, "downloaded: %s | %s" % (ut, uc)
        tried.append("%s   [no valid 13x13 covariance found]" % ut)
    raise RuntimeError(
        "DESI DR2 BAO data (13 points) + full 13x13 covariance could not "
        "be obtained.\nAttempted data sources:\n  " + "\n  ".join(tried)
        + "\nProvide the files explicitly:\n"
          "    python %s <desi_dr2_bao.txt> <desi_dr2_bao.cov>\n"
          "(txt lines: 'z_eff  observable  value [error]'; "
          "cov: full 13x13 matrix)" % sys.argv[0])


# ----------------------------------------------------------------------
# 5. Main
# ----------------------------------------------------------------------

def main():
    line = "=" * 74
    dash = "-" * 74
    print(line)
    print("TASK 1 (BAO, DESI DR2): Delta chi2 = chi2(bimetric) - chi2(LCDM)")
    print(line)
    print("python %s | numpy %s | scipy %s"
          % (sys.version.split()[0], np.__version__, scipy.__version__))
    print("fixed: theta=%.4f deg, m_FP/H0=%.4f, alpha=%.3g, beta=%.3g, "
          "Om_L=%.4f, Om_r=%.5e ; LCDM control Om_m=%.4f"
          % (THETA_DEG, M_FP_OVER_H0, ALPHA, BETA, OM_LAM, OM_R, OM_M_LCDM))

    # ---- bimetric closure ----
    om_m_bim = solve_Om_m_closure()
    y0 = y_physical(om_m_bim + OM_R)
    E2b = make_E2_bim(om_m_bim)
    print(dash)
    print("bimetric closure E^2(0) = 1:")
    print("  Om_m (from closure) = %.10f" % om_m_bim)
    print("  y(0)                = %.12f   [1 - y(0) = %.6e]" % (y0, 1.0 - y0))
    print("  Om_DE(0)            = %.10f" % Om_DE_of_y(y0))
    print("  check: E^2(0) - 1   = %+.3e" % (E2b(0.0) - 1.0))
    print("  LCDM control: Om_L  = 1 - Om_m - Om_r = %.10f"
          % (1.0 - OM_M_LCDM - OM_R))

    # ---- data ----
    print(dash)
    rows, C, src = load_data()
    validate(rows, C)
    n = len(rows)
    zs = np.array([r[0] for r in rows], dtype=float)
    obs = [r[1] for r in rows]
    d = np.array([r[2] for r in rows], dtype=float)
    sig = np.sqrt(np.diag(C))
    print("data source      : %s" % src)
    print("covariance choice: FULL 13x13 from CobayaSampler/bao_data "
          "(desi_bao_dr2)")
    print("parsed data (%d points):" % n)
    for i in range(n):
        note = ""
        if rows[i][3] is not None:
            note = "   [txt err %.5f ; |err-cov|/cov = %.1e]" % (
                rows[i][3], abs(rows[i][3] - sig[i]) / sig[i])
        print("  %2d  z=%8.4f  %s  data=%9.5f  cov_sigma=%8.5f%s"
              % (i + 1, zs[i], obs[i], d[i], sig[i], note))

    # ---- theory vectors and profile minimization over f ----
    v_bim = theory_vector(E2b, zs, obs)
    v_lcdm = theory_vector(E2_lcdm, zs, obs)
    chi2_b, f_b, sf_b = profile_f(d, C, v_bim)
    chi2_l, f_l, sf_l = profile_f(d, C, v_lcdm)
    chk_b = minimize_scalar(lambda f: chi2_at_f(f, d, C, v_bim),
                            bounds=(0.1, 1000.0), method="bounded",
                            options={"xatol": 1e-9})
    chk_l = minimize_scalar(lambda f: chi2_at_f(f, d, C, v_lcdm),
                            bounds=(0.1, 1000.0), method="bounded",
                            options={"xatol": 1e-9})

    print(dash)
    print("predictions at best-fit f (for the log only):")
    print("   #     z_eff  obs      data     bimetric       LCDM")
    for i in range(n):
        print("  %2d  %8.4f  %s  %9.5f  %10.5f  %10.5f"
              % (i + 1, zs[i], obs[i], d[i], f_b * v_bim[i], f_l * v_lcdm[i]))
    print(dash)
    print("f-minimization cross-check (analytic vs scipy.minimize_scalar):")
    print("  bimetric: f = %.8f vs %.8f ; chi2 = %.8f vs %.8f"
          % (f_b, chk_b.x, chi2_b, chk_b.fun))
    print("  LCDM    : f = %.8f vs %.8f ; chi2 = %.8f vs %.8f"
          % (f_l, chk_l.x, chi2_l, chk_l.fun))

    dof = n - 1
    print(line)
    print("RESULTS (TASK 1)")
    print("  Om_m (bimetric, closure E^2(0)=1)  = %.8f" % om_m_bim)
    print("  y(0) (bimetric)                    = %.10f" % y0)
    print("  chi2_bimetric = %.6f   (best f = %.6f +- %.6f)"
          % (chi2_b, f_b, sf_b))
    print("  chi2_LCDM     = %.6f   (best f = %.6f +- %.6f)"
          % (chi2_l, f_l, sf_l))
    print("  Delta chi2 = chi2_bimetric - chi2_LCDM = %+.6f"
          % (chi2_b - chi2_l))
    print("  chi2/dof (dof = %d - 1 = %d): bimetric %.5f , LCDM %.5f"
          % (n, dof, chi2_b / dof, chi2_l / dof))
    print(line)


if __name__ == "__main__":
    main()
