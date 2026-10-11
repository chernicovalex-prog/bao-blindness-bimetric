#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cmb_theta_scan.py
=================
Независимая («слепая») реализация CMB theta-scan для биметрической космологии
по спецификации задания (rev 42). Сжатый CMB-likelihood (Bansal & Huterer 2025):
(R, l_a, omega_b) + гауссова ковариация.

СТАТУС: NOT RUN
    Файл написан, но автором НЕ запускался. Все числовые результаты даёт
    фактический запуск; «ожидаемых чисел» в коде нет.

ВХОДНЫЕ ФАЙЛЫ: НЕ ТРЕБУЮТСЯ
    Все входы — константы ниже, взятые только из текста задания:
    data vector, ковариация (x1e-8), z_rec=1090, omega_gamma=2.469e-5,
    N_eff=3.046, Sum_m_nu=0.06 eV, alpha=0, beta=3, m_FP/H0=25, Omega_L=0.699,
    theta=6.1/15/30 deg, LCDM-контроль (0.309, 0.691),
    границы подгонки h in [0.55,0.88], omega_b in [0.015,0.025].
    Зависимости: python>=3.8, numpy, scipy. Сеть и данные не нужны.

ЗАПУСК
    python3 cmb_theta_scan.py                      # основной прогон + таблица
    python3 cmb_theta_scan.py --check              # только самопроверки
    python3 cmb_theta_scan.py --fast               # быстрее/грубее
    python3 cmb_theta_scan.py --no-grid            # без coarse-grid прескана
    python3 cmb_theta_scan.py --cs-mode literal    # A/B-тест риска R3
    python3 cmb_theta_scan.py --lcdm-mode closure_Omega_L   # A/B-тест R4
        (варианты --lcdm-mode: literal | closure_Omega_L | closure_Omega_m)
    python3 cmb_theta_scan.py --rs-zmax 1e9        # A/B-тест риска R5

МОДЕЛЬ (по спецификации)
    Квартика y(z)=a_f/a_g (коэффициенты по Hogas & Mortsell 2025, Eq. A13;
    m^2 в единицах H0^2):
        a0 = -(1/3) cos^2(th) m^2 (1+2a+be)
        a1 = Om_tot + Om_L + m^2 [ cos^2(th)(a+be) - sin^2(th)(1+a+be/3) ]
        a2 = m^2 [ -cos^2(th) be + sin^2(th)(1+2a+be) ]
        a3 = -[ Om_L + (1/3) m^2 ( cos^2(th)(-1+a-be) + 3 sin^2(th)(a+be) ) ]
        a4 = (1/3) sin^2(th) m^2 be          (a=alpha, be=beta)
    Порядок коэффициентов (риск R1): P(y) = a0 + a1*y + a2*y^2 + a3*y^3 + a4*y^4.
    Обоснование: тождество P(1)=sum(a_i)+Om_tot=Om_tot (проверяется в --check
    на случайных параметрах) и a0<0 => на скобке [1e-6,1] есть смена знака при
    любом z. При обратном порядке (a0*y^4+...+a4) смены знака на [1e-6,1] нет
    ни при каком z, и «Брент на [1e-6,1]» из задания теряет смысл.

    Om_tot(z) = Om_m (1+z)^3 + Om_r (1+z)^4          (массивные nu НЕ входят)
    Om_DE(z)  = Om_L - sin^2(th) m^2 (1-y)[1 + a(1-y) + (be/3)(1-y)^2]
    E^2(z)    = Om_m(1+z)^3 + Om_r(1+z)^4 + Om_DE(z) + rho_nu/rhocrit0
    rho_nu:   (1+z)<=131 -> Om_nu (1+z)^3;  (1+z)>131 -> Om_nu (1+z)^4/131
              (непрерывно на стыке)
    omega_r = omega_gamma(1+(7/8)N_eff(4/11)^(4/3)); omega_nu = 0.06/93.14.

    Самосогласованность: Om_m из E^2(0)=1 (Брент, Om_m in [0.01,0.99]);
    Om_r=omega_r/h^2, Om_nu=omega_nu/h^2 пересчитываются на каждом h,
    Om_m — заново.

    Расстояния (безразмерные, в единицах c/H0 — единственная совместная
    нормировка, при которой R~O(1), l_a~O(300)):
      D_M = int_0^1090 dz/E(z)                  (адаптивная квадратура)
      r_s = int c_s dz/E, по a с якобианом 1/a^2:
            a от 1/(1+z_max) до 1/1091 <=> z от z_max до 1090  (риск R5)
      c_s = [3(1+R_b)]^(-1/2), R_b = 3 omega_b/(4 omega_gamma (1+z))
            (прочтение риска R3; --cs-mode literal — буквальное из задания)
      R   = sqrt(Om_m+Om_nu)*D_M                (префактор по B&H Eq. 4;
                                                  Om_m без массивных nu)
      l_a = pi D_M / r_s

    chi^2 = dV^T C^{-1} dV, C = 1e-8*(матрица из задания), dV = V_model - V_data.
    Подгонка: Nelder-Mead по (h, omega_b) в [0.55,0.88]x[0.015,0.025],
    мультстарт + coarse-grid прескан + полировка.

Ожидаемое время (оценка, НЕ измерено): основной прогон — минуты на обычном
ноутбуке; --fast — быстрее. Если медленно: --fast / --no-grid.
"""

import argparse
import math
import time

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq, minimize

# =============================================================================
# 0. Входные данные (только из текста задания; внешних файлов нет)
# =============================================================================
Z_REC     = 1090.0
A_REC     = 1.0 / 1091.0            # a_rec = 1/1091  <=>  z_rec = 1090
Z_RS_MAX  = 1.0e7                   # верхний предел интеграла r_s (по заданию)

V_DATA = np.array([1.7504, 301.77, 0.022371])            # (R, l_a, omega_b)
COV    = 1.0e-8 * np.array([[1559.83,  -1325.41,   -36.45],
                            [-1325.41, 714691.80,  269.77],
                            [-36.45,    269.77,     2.10]])
COV_INV = np.linalg.inv(COV)
SIG     = np.sqrt(np.diag(COV))     # sigma_R, sigma_la, sigma_wb (для отчёта)

OMEGA_GAMMA  = 2.469e-5
N_EFF        = 3.046
OMEGA_R_BAR  = OMEGA_GAMMA*(1.0 + 0.875*N_EFF*(4.0/11.0)**(4.0/3.0))  # omega_r
OMEGA_NU_BAR = 0.06/93.14                                            # omega_nu

ALPHA       = 0.0
BETA        = 3.0
M2          = 25.0**2                 # (m_FP/H0)^2
OMEGA_L_BIM = 0.699
THETAS_DEG  = (6.1, 15.0, 30.0)

OMEGA_M_LCDM = 0.309
OMEGA_L_LCDM = 0.691

H_LO, H_HI   = 0.55, 0.88
WB_LO, WB_HI = 0.015, 0.025

Y_LO        = 1.0e-6                 # нижний край скобки квартики (по заданию)
NU_W_SWITCH = 131.0                  # 1+z, ниже которой nu нерелятивистские
NU_Z_SWITCH = NU_W_SWITCH - 1.0      # z = 130

SEEDS = [(0.68, 0.0224), (0.60, 0.020), (0.74, 0.0224),
         (0.80, 0.019), (0.62, 0.024), (0.86, 0.0224)]


class Cfg(object):
    """Глобальный конфиг (ставится из аргументов командной строки)."""
    lcdm_mode = "literal"     # literal | closure_Omega_L | closure_Omega_m (R4)
    cs_mode   = "standard"    # standard | literal                           (R3)
    rs_zmax   = Z_RS_MAX      # верхний предел r_s по z                       (R5)


# =============================================================================
# 1. Квартика и модели
# =============================================================================
def quartic_coeffs(s2, c2, alpha, beta, m2, Omega_L):
    """Коэффициенты квартики из задания. Возврат: (a0, a1_off, a2, a3, a4);
    a1 = Omega_tot + a1_off (только a1 зависит от z)."""
    a0 = -(1.0/3.0)*c2*m2*(1.0 + 2.0*alpha + beta)
    a1_off = (Omega_L
              + m2*(c2*(alpha + beta) - s2*(1.0 + alpha + beta/3.0)))
    a2 = m2*(-c2*beta + s2*(1.0 + 2.0*alpha + beta))
    a3 = -(Omega_L + (1.0/3.0)*m2*(c2*(-1.0 + alpha - beta)
                                    + 3.0*s2*(alpha + beta)))
    a4 = (1.0/3.0)*s2*m2*beta
    return a0, a1_off, a2, a3, a4


class BimetricModel(object):
    kind = "bimetric"

    def __init__(self, theta_deg, alpha=ALPHA, beta=BETA, m2=M2,
                 Omega_L=OMEGA_L_BIM):
        self.theta_deg = float(theta_deg)
        th = math.radians(self.theta_deg)
        self.s2 = math.sin(th)**2
        self.c2 = math.cos(th)**2
        self.alpha, self.beta, self.m2 = alpha, beta, m2
        self.Omega_L = Omega_L
        (self.a0, self.a1_off, self.a2, self.a3, self.a4) = quartic_coeffs(
            self.s2, self.c2, alpha, beta, m2, Omega_L)
        self.k_de = self.s2*m2   # множитель при u*(1+a*u+(b/3)u^2), u = 1-y
        # Достаточное условие единственности корня на [0,1]:
        # P'(y) = a1 + 2 a2 y + 3 a3 y^2 + 4 a4 y^3, причём
        # a1 = Om_tot + a1_off >= 0.01 + a1_off при любом допустимом
        # Om_m (>= 0.01 из скобки замыкания) и любом z.
        ys = np.linspace(0.0, 1.0, 2001)
        self.p_prime_min = float(np.min(
            (0.01 + self.a1_off) + 2.0*self.a2*ys
            + 3.0*self.a3*ys**2 + 4.0*self.a4*ys**3))
        self.n_clamped = 0
        self._uni = {}   # кэш Universe по h
        self._dM = {}    # кэш D_M по h (D_M не зависит от omega_b)


class LcdmModel(object):
    kind = "lcdm"

    def __init__(self):
        self.Omega_m_in = OMEGA_M_LCDM
        self.Omega_L_in = OMEGA_L_LCDM
        self.n_clamped = 0
        self._uni = {}
        self._dM = {}


def quartic_root(model, Omega_tot):
    """Физическая ветвь y(z): корень P(y)=0 на [Y_LO, 1] методом Брента.

    P(y) = a0 + a1 y + a2 y^2 + a3 y^3 + a4 y^4, a1 = Om_tot + a1_off.
    Если P(1e-6) >= 0, физический корень опустился ниже 1e-6 (большие z);
    тогда возвращаем 1e-6 (риск R2): (1-y) отличается от 1 не более чем на
    1e-6, а вклад этого члена в E^2 при z >~ 1e3 ничтожен. Число клампов
    накапливается в model.n_clamped и печатается.
    """
    if model.kind != "bimetric":
        raise TypeError("quartic_root: только для биметрических моделей")
    a0, a2, a3, a4 = model.a0, model.a2, model.a3, model.a4
    a1 = Omega_tot + model.a1_off

    def P(y):
        return a0 + y*(a1 + y*(a2 + y*(a3 + y*a4)))

    if P(Y_LO) >= 0.0:
        model.n_clamped += 1
        return Y_LO
    return brentq(P, Y_LO, 1.0, xtol=1.0e-13, maxiter=200)


# =============================================================================
# 2. Вселенная при данном h (замыкание Om_m из E^2(0)=1)
# =============================================================================
class Universe(object):
    def __init__(self, model, h):
        self.model = model
        self.h = float(h)
        self.Omega_r = OMEGA_R_BAR/self.h**2
        self.Omega_nu = OMEGA_NU_BAR/self.h**2
        if model.kind == "bimetric":
            self.Omega_L = model.Omega_L
            self.Omega_m = self._solve_Omega_m()
            self.z_clamp = self._find_z_clamp()
        elif Cfg.lcdm_mode == "literal":
            # Буквально из задания: оба параметра фиксированы, тогда
            # E^2(0) = 1 + Om_r + Om_nu (значение печатается в отчёте).
            self.Omega_m = model.Omega_m_in
            self.Omega_L = model.Omega_L_in
            self.z_clamp = None
        elif Cfg.lcdm_mode == "closure_Omega_L":
            self.Omega_L = model.Omega_L_in
            self.Omega_m = self._solve_Omega_m()
            self.z_clamp = None
        elif Cfg.lcdm_mode == "closure_Omega_m":
            self.Omega_m = model.Omega_m_in
            self.Omega_L = 1.0 - self.Omega_m - self.Omega_r - self.Omega_nu
            self.z_clamp = None
        else:
            raise ValueError("unknown lcdm_mode: %r" % Cfg.lcdm_mode)

    def _DE_at_z0(self, Omega_m):
        """Omega_DE(z=0) при пробном Omega_m (для условия E^2(0)=1)."""
        m = self.model
        if m.kind == "lcdm":
            return self.Omega_L
        y0 = quartic_root(m, Omega_m + self.Omega_r)  # Om_tot(0)=Om_m+Om_r
        u = 1.0 - y0
        return m.Omega_L - m.k_de*u*(1.0 + m.alpha*u + (m.beta/3.0)*u*u)

    def _solve_Omega_m(self):
        """Omega_m из E^2(0)=1; метод Брента, Omega_m in [0.01, 0.99]."""
        def g(Om):
            return (Om + self.Omega_r + self.Omega_nu
                    + self._DE_at_z0(Om) - 1.0)
        return brentq(g, 0.01, 0.99, xtol=1.0e-13, maxiter=200)

    def _find_z_clamp(self):
        """z, в котором корень y достигает 1e-6 (ниже — кламп R2)."""
        m = self.model

        def hfun(z):
            w = 1.0 + z
            Ot = self.Omega_m*w**3 + self.Omega_r*w**4
            a1 = Ot + m.a1_off
            return m.a0 + Y_LO*(a1 + Y_LO*(m.a2 + Y_LO*(m.a3 + Y_LO*m.a4)))

        try:
            return brentq(hfun, 10.0, Z_RS_MAX, xtol=1.0e-6, maxiter=200)
        except ValueError:
            return None


def get_universe(model, h):
    if h in model._uni:
        return model._uni[h]
    uni = Universe(model, h)
    model._uni[h] = uni
    return uni


# =============================================================================
# 3. E^2(z), c_s, расстояния
# =============================================================================
def Omega_DE(model, y):
    u = 1.0 - y
    return (model.Omega_L
            - model.k_de*u*(1.0 + model.alpha*u + (model.beta/3.0)*u*u))


def nu_density(uni, w):
    if w <= NU_W_SWITCH:
        return uni.Omega_nu*w**3
    return uni.Omega_nu*w**4/NU_W_SWITCH


def E2(uni, z):
    w = 1.0 + z
    Om_tot = uni.Omega_m*w**3 + uni.Omega_r*w**4
    if uni.model.kind == "bimetric":
        DE = Omega_DE(uni.model, quartic_root(uni.model, Om_tot))
    else:
        DE = uni.Omega_L
    return Om_tot + DE + nu_density(uni, w)


def c_sound(z, omega_b):
    """c_s в единицах c. Риск R3: два прочтения формулы задания."""
    w = 1.0 + z
    R_b = 3.0*omega_b/(4.0*OMEGA_GAMMA*w)
    if Cfg.cs_mode == "standard":
        return 1.0/math.sqrt(3.0*(1.0 + R_b))
    return 1.0/math.sqrt(3.0*(1.0 + R_b/w))     # 'literal'


def _pts_inside(pts, lo, hi):
    return sorted(float(p) for p in pts if lo < float(p) < hi)


def comoving_D(uni):
    """D_M = int_0^1090 dz/E(z), адаптивная квадратура; единицы c/H0.
    Разбиение в точке переключения nu (z=130) и в точке клампа (если есть)."""
    m = uni.model
    if uni.h in m._dM:
        return m._dM[uni.h]

    def f(z):
        return 1.0/math.sqrt(E2(uni, z))

    low = _pts_inside([1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0],
                      0.0, NU_Z_SWITCH)
    I = quad(f, 0.0, NU_Z_SWITCH, points=low, limit=400,
             epsabs=1e-10, epsrel=1e-10)[0]
    hp = [200.0, 400.0, 700.0, 1000.0]
    if uni.z_clamp is not None:
        hp.append(uni.z_clamp)
    hp = _pts_inside(hp, NU_Z_SWITCH, Z_REC)
    I += quad(f, NU_Z_SWITCH, Z_REC, points=hp, limit=400,
              epsabs=1e-10, epsrel=1e-10)[0]
    m._dM[uni.h] = I
    return I


def sound_horizon(uni, omega_b):
    """r_s = int c_s dz/E; по переменной a с якобианом 1/a^2:
    a от 1/(1+rs_zmax) до 1/1091  <=>  z от rs_zmax до 1090. Единицы c/H0."""
    a_lo = 1.0/(1.0 + Cfg.rs_zmax)

    def integrand(a):
        z = 1.0/a - 1.0
        return c_sound(z, omega_b)/(a*a*math.sqrt(E2(uni, z)))

    pts = [a_lo*9.0, 1e-6, 1e-5, 1e-4, 2e-4, 5e-4]
    if uni.z_clamp is not None:
        pts.append(1.0/(1.0 + uni.z_clamp))
    pts = _pts_inside(pts, a_lo, A_REC)
    return quad(integrand, a_lo, A_REC, points=pts, limit=400,
                epsabs=1e-12, epsrel=1e-10)[0]


def predict(model, h, omega_b):
    uni = get_universe(model, h)
    DM = comoving_D(uni)
    rs = sound_horizon(uni, omega_b)
    R = math.sqrt(uni.Omega_m + uni.Omega_nu)*DM
    l_a = math.pi*DM/rs
    return R, l_a


# =============================================================================
# 4. chi^2 и подгонка (Nelder-Mead, мультстарт)
# =============================================================================
def make_chi2(model):
    state = {"n": 0, "fail": 0}

    def chi2(x):
        h = float(np.clip(x[0], H_LO, H_HI))
        wb = float(np.clip(x[1], WB_LO, WB_HI))
        state["n"] += 1
        try:
            R, la = predict(model, h, wb)
            d = np.array([R, la, wb]) - V_DATA
            return float(d @ COV_INV @ d)
        except Exception as exc:
            state["fail"] += 1
            print("  [warn] сбой вычисления при h=%.4f, wb=%.5f: %r"
                  % (h, wb, exc))
            return 1.0e10

    chi2.state = state
    return chi2


def _nm(func, x0, bounds, opts):
    try:
        return minimize(func, x0=x0, method="Nelder-Mead", bounds=bounds,
                        options=opts)
    except TypeError:                      # старый scipy: у NM нет bounds
        return minimize(func, x0=x0, method="Nelder-Mead", options=opts)


def fit(model, fast=False, grid_prescan=True):
    chi2 = make_chi2(model)
    seeds = list(SEEDS[:3] if fast else SEEDS)
    if grid_prescan:
        hv = np.linspace(H_LO + 0.03, H_HI - 0.03, 7)
        wv = np.linspace(WB_LO + 0.001, WB_HI - 0.001, 5)
        bestg = None
        for h in hv:
            for wb in wv:
                v = chi2((float(h), float(wb)))
                if bestg is None or v < bestg[0]:
                    bestg = (v, (float(h), float(wb)))
        if bestg is not None:
            seeds.append(bestg[1])
    bounds = [(H_LO, H_HI), (WB_LO, WB_HI)]
    opts = dict(xatol=(1e-3 if fast else 1e-4),
                fatol=(1e-4 if fast else 1e-6),
                maxiter=2000, maxfev=4000)
    best = None
    for s in seeds:
        res = _nm(chi2, np.array(s, dtype=float), bounds, opts)
        if best is None or res.fun < best.fun:
            best = res
    if best is None:
        raise RuntimeError("fit: все вычисления chi2 завершились сбоем")
    for _ in range(2):                     # полировка
        res = _nm(chi2, np.asarray(best.x, dtype=float), bounds,
                  dict(xatol=1e-5, fatol=1e-8, maxiter=2000, maxfev=4000))
        if res.fun < best.fun:
            best = res
    h = float(np.clip(best.x[0], H_LO, H_HI))
    wb = float(np.clip(best.x[1], WB_LO, WB_HI))
    uni = get_universe(model, h)
    R, la = predict(model, h, wb)
    tag = "LCDM" if model.kind == "lcdm" else "%g deg" % model.theta_deg
    return dict(
        model=model, tag=tag, chi2=float(best.fun), h=h, wb=wb,
        Omega_m=uni.Omega_m, Omega_L=uni.Omega_L,
        R=R, la=la, DM=comoving_D(uni), rs=sound_horizon(uni, wb),
        E2_0=E2(uni, 0.0),
        dR_sig=(R - V_DATA[0])/SIG[0],
        dla_sig=(la - V_DATA[1])/SIG[1],
        dwb_sig=(wb - V_DATA[2])/SIG[2],
        n_clamped=model.n_clamped, n_eval=chi2.state["n"],
        n_fail=chi2.state["fail"],
        on_bound=(min(h - H_LO, H_HI - h) < 1e-3
                  or min(wb - WB_LO, WB_HI - wb) < 1e-5),
    )


# =============================================================================
# 5. Самопроверки (--check)
# =============================================================================
def _gl_panels(panels, n=64):
    """Композитная квадратура Гаусса-Лежандра по списку краёв панелей."""
    x, w = np.polynomial.legendre.leggauss(n)
    xs, ws = [], []
    for lo, hi in zip(panels[:-1], panels[1:]):
        xm, xr = 0.5*(lo + hi), 0.5*(hi - lo)
        xs.append(xm + xr*x)
        ws.append(xr*w)
    return np.concatenate(xs), np.concatenate(ws)


def do_checks():
    print("=" * 78)
    print("САМОПРОВЕРКИ (--check); lcdm-mode=%s, cs-mode=%s"
          % (Cfg.lcdm_mode, Cfg.cs_mode))
    print("=" * 78)

    # [1] тождество коэффициентов: P(1) = Omega_tot при любых параметрах
    rng = np.random.default_rng(20250501)
    worst = 0.0
    for _ in range(500):
        th = rng.uniform(0.0, 0.5*math.pi)
        s2, c2 = math.sin(th)**2, math.cos(th)**2
        al = rng.uniform(-1.0, 1.0)
        be = rng.uniform(-1.0, 4.0)
        m2 = 10.0 ** rng.uniform(0.0, 4.0)
        OL = rng.uniform(0.0, 1.0)
        Ot = 10.0 ** rng.uniform(-3.0, 8.0)
        a0, a1o, a2, a3, a4 = quartic_coeffs(s2, c2, al, be, m2, OL)
        S = a0 + (Ot + a1o) + a2 + a3 + a4
        worst = max(worst, abs(S - Ot)/max(1.0, abs(Ot)))
    print("[1] P(1)=Om_tot: max rel dev = %.2e  ->  %s"
          % (worst, "OK" if worst < 1e-10 else "FAIL"))

    # [2] монотонность P на [0,1] => корень на скобке единственен
    for td in THETAS_DEG:
        m = BimetricModel(td)
        print("[2] theta=%5.1f deg: min_y P'(y) (худший случай a1) = %+12.3f"
              "  ->  %s" % (td, m.p_prime_min,
                            "OK (корень единственен)" if m.p_prime_min > 0.0
                            else "ВНИМАНИЕ: возможны несколько корней"))

    # [3] Brent против numpy.roots (в зоне, где корень > 1e-6)
    for td in THETAS_DEG:
        m = BimetricModel(td)
        uni = Universe(m, 0.68)
        worst, bad = 0.0, 0
        for z in (0.0, 0.1, 0.5, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0):
            w = 1.0 + z
            Ot = uni.Omega_m*w**3 + uni.Omega_r*w**4
            a1 = Ot + m.a1_off
            Plo = m.a0 + Y_LO*(a1 + Y_LO*(m.a2 + Y_LO*(m.a3 + Y_LO*m.a4)))
            if not (Plo < 0.0):
                bad += 1
                continue
            yb = quartic_root(m, Ot)
            roots = np.roots([m.a4, m.a3, m.a2, a1, m.a0])
            inside = sorted(r.real for r in roots
                            if abs(r.imag) <= 1e-8*max(1.0, abs(r))
                            and 1e-12 < r.real < 1.0)
            if len(inside) != 1:
                bad += 1
            else:
                worst = max(worst, abs(inside[0] - yb))
        print("[3] theta=%5.1f deg: |brent - numpy.roots| = %.2e, аномалий = %d"
              "  ->  %s" % (td, worst, bad,
                            "OK" if bad == 0 and worst < 1e-7 else "FAIL"))

    # [4] непрерывность E^2 на стыке 1+z=131
    m = BimetricModel(15.0)
    uni = Universe(m, 0.68)
    e = 1e-9
    lo, hi = E2(uni, NU_Z_SWITCH - e), E2(uni, NU_Z_SWITCH + e)
    rel = abs(lo - hi)/hi
    print("[4] непрерывность E^2 на 1+z=131: rel jump = %.2e  ->  %s"
          % (rel, "OK" if rel < 1e-8 else "FAIL"))

    # [5] предел theta -> 0: биметрия вырождается в LCDM с Om_L = 0.699
    m0 = BimetricModel(0.0)
    u0 = Universe(m0, 0.68)
    dev_m = abs(u0.Omega_m - (1.0 - u0.Omega_r - u0.Omega_nu - 0.699))
    worst_e = 0.0
    for z in (0.0, 1.0, 10.0, 130.0, 500.0, Z_REC, 1e4, 1e6):
        w = 1.0 + z
        nu = (u0.Omega_nu*w**3 if w <= NU_W_SWITCH
              else u0.Omega_nu*w**4/NU_W_SWITCH)
        direct = u0.Omega_m*w**3 + u0.Omega_r*w**4 + 0.699 + nu
        worst_e = max(worst_e, abs(E2(u0, z) - direct)/max(1.0, direct))
    print("[5] theta=0: |Om_m - (1-Om_r-Om_nu-0.699)| = %.2e, max rel dev E^2 ="
          " %.2e  ->  %s" % (dev_m, worst_e,
                             "OK" if dev_m < 1e-9 and worst_e < 1e-12
                             else "FAIL"))

    # [6] квадратуры: scipy.quad против композитной Гаусса-Лежандра
    wb0 = 0.022371
    mm = BimetricModel(15.0)
    uu = Universe(mm, 0.68)
    DM_q = comoving_D(uu)
    pan = [0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 20.0, 30.0,
           50.0, 80.0, 120.0, NU_Z_SWITCH, 200.0, 400.0, 700.0, 1000.0, Z_REC]
    xs, ws = _gl_panels(pan, 64)
    DM_g = float(np.dot(ws, [1.0/math.sqrt(E2(uu, float(z))) for z in xs]))
    rs_q = sound_horizon(uu, wb0)
    pana = [1.0/(1.0 + Cfg.rs_zmax), 1e-6, 1e-5, 1e-4, 2e-4, 4e-4, A_REC]
    xa, wa_ = _gl_panels(pana, 64)
    rs_g = float(np.dot(wa_, [c_sound(1.0/float(a) - 1.0, wb0)
                              /(a*a*math.sqrt(E2(uu, 1.0/float(a) - 1.0)))
                              for a in xa]))
    print("[6] D_M: quad=%.10f GL=%.10f (rel %.2e); r_s: quad=%.10f GL=%.10f"
          " (rel %.2e)  ->  %s"
          % (DM_q, DM_g, abs(DM_q - DM_g)/DM_q,
             rs_q, rs_g, abs(rs_q - rs_g)/rs_g,
             "OK" if abs(DM_q - DM_g)/DM_q < 1e-8
             and abs(rs_q - rs_g)/rs_q < 1e-8 else "FAIL"))

    # [7] ковариация
    corr = COV/np.outer(SIG, SIG)
    print("[7] sigma = (%.6f, %.6f, %.2e); rho(R,la)=%.4f rho(R,wb)=%.4f "
          "rho(la,wb)=%.4f" % (SIG[0], SIG[1], SIG[2],
                               corr[0, 1], corr[0, 2], corr[1, 2]))
    try:
        np.linalg.cholesky(COV)
        print("    ковариация положительно определена: OK")
    except np.linalg.LinAlgError:
        print("    ВНИМАНИЕ: ковариация не PD")

    # [8] зона клампа (риск R2): где корень y опускается ниже 1e-6
    for td in THETAS_DEG:
        mm = BimetricModel(td)
        uu = Universe(mm, 0.68)
        zc = uu.z_clamp
        print("[8] theta=%5.1f deg: корень y достигает 1e-6 при z = %s"
              % (td, ("%.1f" % zc) if zc is not None else "нет в [10, 1e7]"))
    print("=" * 78)


# =============================================================================
# 6. Основной прогон
# =============================================================================
def run(fast=False, grid=True):
    print("data vector V =", V_DATA, "  sigma =", SIG)
    models = [BimetricModel(t) for t in THETAS_DEG] + [LcdmModel()]
    results = []
    t00 = time.time()
    for m in models:
        t0 = time.time()
        print("\n" + "-" * 78)
        if m.kind == "bimetric":
            print("Модель: theta = %g deg (alpha=%g, beta=%g, m_FP/H0=%g,"
                  " Omega_L=%g)" % (m.theta_deg, m.alpha, m.beta,
                                    math.sqrt(m.m2), m.Omega_L))
            print("  гарантия единственности корня: min P' = %+.3f"
                  % m.p_prime_min)
        else:
            print("Модель: LCDM-контроль (Omega_m=%g, Omega_L=%g, режим '%s')"
                  % (m.Omega_m_in, m.Omega_L_in, Cfg.lcdm_mode))
        smoke = make_chi2(m)
        print("  smoke chi2(h=0.68, omega_b=0.022371) = %.4f"
              % smoke((0.68, 0.022371)))
        info = fit(m, fast=fast, grid_prescan=grid)
        info["elapsed"] = time.time() - t0
        results.append(info)
        print("  chi2    = %.6f" % info["chi2"])
        print("  h_best  = %.5f   omega_b = %.6f" % (info["h"], info["wb"]))
        print("  Omega_m = %.6f   (Omega_L = %.6f, E^2(0) = %.10f)"
              % (info["Omega_m"], info["Omega_L"], info["E2_0"]))
        print("  R = %.5f, l_a = %.4f  (отклонения: %+.2f, %+.2f, %+.2f сигм)"
              % (info["R"], info["la"], info["dR_sig"], info["dla_sig"],
                 info["dwb_sig"]))
        print("  D_M = %.6f, r_s = %.8f" % (info["DM"], info["rs"]))
        print("  вызовов с клампом y<1e-6: %d; вычислений chi2: %d (сбоев %d);"
              " время %.1f c" % (info["n_clamped"], info["n_eval"],
                                 info["n_fail"], info["elapsed"]))
        if info["on_bound"]:
            print("  ВНИМАНИЕ: оптимум на границе диапазона (h или omega_b).")
    print("\n" + "=" * 78)
    print("ИТОГОВАЯ ТАБЛИЦА:  theta | chi2 | h_best | Omega_m")
    print("=" * 78)
    print("%-12s %14s %12s %12s" % ("theta", "chi2", "h_best", "Omega_m"))
    for r in results:
        print("%-12s %14.6f %12.5f %12.5f"
              % (r["tag"], r["chi2"], r["h"], r["Omega_m"]))
    print("конфигурация: lcdm-mode=%s, cs-mode=%s, rs-zmax=%g, fast=%s, "
          "grid=%s" % (Cfg.lcdm_mode, Cfg.cs_mode, Cfg.rs_zmax, fast, grid))
    print("полное время: %.1f c" % (time.time() - t00))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="CMB theta-scan для биметрической космологии "
                    "(подробности — в docstring файла)")
    ap.add_argument("--fast", action="store_true",
                    help="урезанный мультстарт и грубые допуски NM")
    ap.add_argument("--no-grid", action="store_true",
                    help="отключить coarse-grid прескан")
    ap.add_argument("--check", action="store_true",
                    help="только самопроверки, без фитов")
    ap.add_argument("--lcdm-mode", default="literal",
                    choices=["literal", "closure_Omega_L", "closure_Omega_m"],
                    help="режим LCDM-контроля (риск R4)")
    ap.add_argument("--cs-mode", default="standard",
                    choices=["standard", "literal"],
                    help="прочтение формулы c_s (риск R3)")
    ap.add_argument("--rs-zmax", type=float, default=Z_RS_MAX,
                    help="верхний предел r_s по z (риск R5)")
    args = ap.parse_args(argv)
    Cfg.lcdm_mode = args.lcdm_mode
    Cfg.cs_mode = args.cs_mode
    Cfg.rs_zmax = args.rs_zmax

    print("=" * 78)
    print("CMB theta-scan (биметрия), независимая реализация по спецификации.")
    print("СТАТУС КОДА: доставлен с пометкой NOT RUN; все числа ниже")
    print("генерируются ТЕКУЩИМ запуском.")
    print("numpy %s, scipy %s; lcdm-mode=%s, cs-mode=%s, rs-zmax=%g"
          % (np.__version__, scipy.__version__,
             Cfg.lcdm_mode, Cfg.cs_mode, Cfg.rs_zmax))
    print("=" * 78)

    if args.check:
        do_checks()
        return 0
    run(fast=args.fast, grid=not args.no_grid)
    return 0


if __name__ == "__main__":
    main()
