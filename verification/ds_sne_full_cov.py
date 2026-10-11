#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
独立复现 Pantheon+ Δχ² 分析：双度量宇宙学背景 vs ΛCDM。
数据文件需放在当前目录：
  Pantheon+SH0ES.dat
  Pantheon+SH0ES_STAT+SYS.cov
"""

import numpy as np
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
import sys
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ----------------------------------------------------------------------
# 1. 读取数据
# ----------------------------------------------------------------------
dat_file = "Pantheon+SH0ES.dat"
cov_file = "Pantheon+SH0ES_STAT+SYS.cov"

# 列索引（0-based）：
# 0 CID, 1 IDSURVEY, 2 zHD, 3 zHDERR, 4 zCMB, 5 zCMBERR,
# 6 zHEL, 7 zHELERR, 8 m_b_corr, 9 m_b_corr_err_DIAG
try:
    data = np.loadtxt(dat_file, comments="#", usecols=(4, 6, 8, 9))
except Exception as e:
    print(f"ERROR: 无法读取 {dat_file}: {e}")
    sys.exit(1)

zCMB = data[:, 0]
zHEL = data[:, 1]
mb   = data[:, 2]
mb_err = data[:, 3]

# 筛选 zCMB > 0.01
mask = zCMB > 0.01
zCMB = zCMB[mask]
zHEL = zHEL[mask]
mb   = mb[mask]
mb_err = mb_err[mask]
N_data = len(zCMB)
print(f"筛选后数据点数 N = {N_data}")

# 读取全协方差
try:
    cov_flat = np.loadtxt(cov_file)
    N_cov = int(cov_flat[0])
    if len(cov_flat) != 1 + N_cov * N_cov:
        raise ValueError(f"协方差文件长度异常: {len(cov_flat)} != 1 + {N_cov}^2")
    C_full = cov_flat[1:].reshape(N_cov, N_cov)
except Exception as e:
    print(f"ERROR: 无法读取 {cov_file}: {e}")
    sys.exit(1)

if N_cov != len(mask):
    print(f"警告: 协方差维度 {N_cov} 与数据行数 {len(mask)} 不一致")
# 应用掩码到协方差
C_full = C_full[np.ix_(mask, mask)]
C_diag = np.diag(mb_err**2)

# ----------------------------------------------------------------------
# 2. 模型常数
# ----------------------------------------------------------------------
c_over_H0 = 1.0          # 任意常数，边缘化 M 后不影响 χ²
Omega_r = 9.08e-5
Omega_Lambda = 0.699
theta_deg = 6.1
theta = np.deg2rad(theta_deg)
m_fp_over_H0 = 25.0
alpha = 0.0
beta = 3.0

m2 = m_fp_over_H0**2
cos2 = np.cos(theta)**2
sin2 = np.sin(theta)**2

# ----------------------------------------------------------------------
# 3. ΛCDM 模型
# ----------------------------------------------------------------------
def E2_lcdm(z, Omega_m):
    Omega_L = 1.0 - Omega_m - Omega_r
    return Omega_m * (1.0 + z)**3 + Omega_r * (1.0 + z)**4 + Omega_L

def E_lcdm(z, Omega_m):
    e2 = E2_lcdm(z, Omega_m)
    if e2 <= 0:
        return np.nan
    return np.sqrt(e2)

# ----------------------------------------------------------------------
# 4. 双度量模型
# ----------------------------------------------------------------------
def get_y_bimetric(z, Omega_m):
    """解四次方程，返回物理根 y ∈ (0,1]；若有多个，取最大。"""
    Omega_tot = Omega_m * (1.0 + z)**3 + Omega_r * (1.0 + z)**4
    a0 = -(1.0/3.0) * cos2 * m2 * (1.0 + 2.0*alpha + beta)
    a1 = Omega_tot + Omega_Lambda + m2 * (cos2*(alpha + beta) - sin2*(1.0 + alpha + beta/3.0))
    a2 = m2 * (-cos2*beta + sin2*(1.0 + 2.0*alpha + beta))
    a3 = -(Omega_Lambda + (1.0/3.0)*m2*(cos2*(-1.0 + alpha - beta) + 3.0*sin2*(alpha + beta)))
    a4 = (1.0/3.0) * sin2 * m2 * beta
    coeffs = [a4, a3, a2, a1, a0]
    roots = np.roots(coeffs)
    real_roots = roots[np.isreal(roots)].real
    valid = real_roots[(real_roots > 0.0) & (real_roots <= 1.0)]
    if len(valid) == 0:
        # 放宽一点边界，防止数值误差
        valid = real_roots[(real_roots > 0.0) & (real_roots <= 1.0 + 1e-8)]
    if len(valid) == 0:
        return np.nan
    return np.max(valid)   # 物理根通常最接近 1

def Omega_DE_bimetric(z, y):
    return Omega_Lambda - sin2 * m2 * (1.0 - y) * (1.0 + alpha*(1.0 - y) + (beta/3.0)*(1.0 - y)**2)

def E2_bimetric(z, Omega_m):
    y = get_y_bimetric(z, Omega_m)
    if not np.isfinite(y):
        return np.nan
    return Omega_m * (1.0 + z)**3 + Omega_r * (1.0 + z)**4 + Omega_DE_bimetric(z, y)

def E_bimetric(z, Omega_m):
    e2 = E2_bimetric(z, Omega_m)
    if e2 <= 0 or not np.isfinite(e2):
        return np.nan
    return np.sqrt(e2)

def closure_residual(Omega_m):
    """E²(0) - 1，用于求自洽 Ω_m。"""
    e2_0 = E2_bimetric(0.0, Omega_m)
    if not np.isfinite(e2_0):
        return np.nan
    return e2_0 - 1.0

# 求自洽 Ω_m
try:
    Omega_m_bim = brentq(closure_residual, 0.01, 0.99, xtol=1e-10, rtol=1e-10, maxiter=200)
except Exception as e:
    print(f"ERROR: 双度量闭合求根失败: {e}")
    Omega_m_bim = np.nan

if np.isfinite(Omega_m_bim):
    y0_bim = get_y_bimetric(0.0, Omega_m_bim)
    print(f"双度量自洽 Ω_m = {Omega_m_bim:.10f}")
    print(f"双度量 y(0) = {y0_bim:.10f}")
else:
    y0_bim = np.nan
    print("双度量自洽 Ω_m 求解失败")

# ----------------------------------------------------------------------
# 5. 距离模数
# ----------------------------------------------------------------------
def comoving_distance(z, E_func, Omega_m=None):
    """∫0^z dz'/E(z')，E_func 为 E(z) 函数。"""
    def integrand(zp):
        if Omega_m is None:
            e = E_func(zp)
        else:
            e = E_func(zp, Omega_m)
        if not np.isfinite(e) or e <= 0:
            return np.nan
        return 1.0 / e
    try:
        val, err = quad(integrand, 0.0, z, epsabs=1e-8, epsrel=1e-8, limit=200)
        return val
    except Exception:
        return np.nan

def mu_model(zCMB_arr, zHEL_arr, E_func, Omega_m=None):
    """计算 μ = 5 log10[(1+z_HEL) * (c/H0) * D_C] + 25。"""
    mu = np.full_like(zCMB_arr, np.nan, dtype=float)
    for i, (zc, zh) in enumerate(zip(zCMB_arr, zHEL_arr)):
        dc = comoving_distance(zc, E_func, Omega_m)
        if not np.isfinite(dc):
            mu[i] = np.nan
            continue
        mu[i] = 5.0 * np.log10((1.0 + zh) * c_over_H0 * dc) + 25.0
    return mu

# ----------------------------------------------------------------------
# 6. χ² 边缘化 M
# ----------------------------------------------------------------------
def chi2_marginalized_diag(r, sigma):
    """对角协方差，解析边缘化 M。"""
    w = 1.0 / sigma**2
    sw = np.sum(w)
    swr = np.sum(w * r)
    swr2 = np.sum(w * r**2)
    return swr2 - (swr**2) / sw

def chi2_marginalized_full(r, cho, C_inv_1):
    """全协方差，使用 Cholesky 分解。"""
    C_inv_r = cho_solve(cho, r)
    ones = np.ones_like(r)
    # 注意：C_inv_1 已预先算好 = C^{-1} 1
    num = (ones @ C_inv_r)**2
    den = ones @ C_inv_1
    return r @ C_inv_r - num / den

# 预计算全协方差的 Cholesky 和 C_inv_1
try:
    cho_full = cho_factor(C_full, lower=True, check_finite=False)
    C_inv_1_full = cho_solve(cho_full, np.ones(N_data))
except Exception as e:
    print(f"ERROR: 全协方差 Cholesky 分解失败: {e}")
    cho_full = None
    C_inv_1_full = None

# ----------------------------------------------------------------------
# 7. ΛCDM 控制：网格/最小化
# ----------------------------------------------------------------------
def chi2_lcdm(Omega_m, cov_type):
    mu = mu_model(zCMB, zHEL, E_lcdm, Omega_m)
    if np.any(~np.isfinite(mu)):
        return 1e300
    r = mb - mu
    if cov_type == "diag":
        return chi2_marginalized_diag(r, mb_err)
    else:
        if cho_full is None:
            return 1e300
        return chi2_marginalized_full(r, cho_full, C_inv_1_full)

results = {}

for cov_type in ["diag", "full"]:
    print(f"\n===== 协方差类型: {cov_type} =====")
    # ΛCDM 最佳 Ω_m
    res = minimize_scalar(lambda om: chi2_lcdm(om, cov_type),
                          bounds=(0.05, 0.60), method="bounded",
                          options={"xatol": 1e-6})
    Omega_m_best = res.x
    chi2_min = res.fun
    print(f"ΛCDM 最佳 Ω_m = {Omega_m_best:.10f}")
    print(f"ΛCDM χ²_min = {chi2_min:.6f}")

    # 双度量 χ²
    if np.isfinite(Omega_m_bim):
        mu_bim = mu_model(zCMB, zHEL, E_bimetric, Omega_m_bim)
        if np.any(~np.isfinite(mu_bim)):
            chi2_bim = np.nan
            print("警告: 双度量 μ 存在非有限值")
        else:
            r_bim = mb - mu_bim
            if cov_type == "diag":
                chi2_bim = chi2_marginalized_diag(r_bim, mb_err)
            else:
                if cho_full is None:
                    chi2_bim = np.nan
                else:
                    chi2_bim = chi2_marginalized_full(r_bim, cho_full, C_inv_1_full)
    else:
        chi2_bim = np.nan

    delta_chi2 = chi2_bim - chi2_min if np.isfinite(chi2_bim) else np.nan

    # 自由度：边缘化 M 相当于拟合了一个常数。
    # ΛCDM 还拟合了 Ω_m，所以 n_params = 2；双度量无其他自由参数，n_params = 1。
    dof_lcdm = N_data - 2
    dof_bim  = N_data - 1
    chi2_dof_lcdm = chi2_min / dof_lcdm if dof_lcdm > 0 else np.nan
    chi2_dof_bim  = chi2_bim / dof_bim if (np.isfinite(chi2_bim) and dof_bim > 0) else np.nan

    print(f"双度量 χ² = {chi2_bim:.6f}")
    print(f"Δχ² = χ²_bim - χ²_LCDM = {delta_chi2:.6f}")
    print(f"χ²/dof (ΛCDM) = {chi2_dof_lcdm:.6f}")
    print(f"χ²/dof (bimetric) = {chi2_dof_bim:.6f}")

    results[cov_type] = {
        "Omega_m_best_lcdm": Omega_m_best,
        "chi2_min_lcdm": chi2_min,
        "chi2_bim": chi2_bim,
        "delta_chi2": delta_chi2,
        "chi2_dof_lcdm": chi2_dof_lcdm,
        "chi2_dof_bim": chi2_dof_bim,
    }

# ----------------------------------------------------------------------
# 8. 汇总输出
# ----------------------------------------------------------------------
print("\n================ 汇总 ================")
print(f"双度量自洽 Ω_m = {Omega_m_bim:.10f}")
print(f"双度量 y(0) = {y0_bim:.10f}")
for cov_type in ["diag", "full"]:
    r = results[cov_type]
    print(f"\n--- {cov_type} 协方差 ---")
    print(f"Ω_m^best(ΛCDM) = {r['Omega_m_best_lcdm']:.10f}")
    print(f"χ²_min(ΛCDM)   = {r['chi2_min_lcdm']:.6f}")
    print(f"χ²_bim         = {r['chi2_bim']:.6f}")
    print(f"Δχ²            = {r['delta_chi2']:.6f}")
    print(f"χ²/dof (ΛCDM)  = {r['chi2_dof_lcdm']:.6f}")
    print(f"χ²/dof (bim)   = {r['chi2_dof_bim']:.6f}")
