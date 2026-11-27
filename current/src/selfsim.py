"""Self-similar reduction of the scalar-instability boundary (paper Sec. 6.4).

The quartic divided by s^2 m_FP^2 (Omega_r = 0) is free of m_FP except for
the explicit epsilon = Omega_Lambda/(s^2 m_FP^2) terms (B-parameters are
linear in Omega_Lambda; see the paper, Lemma 1):

    P(y; u, eps) with u = Omega_m(1+z)^3/(s^2 m_FP^2):
    a4 = t*beta
    a3 = -3 t (alpha+beta) + (1-alpha+beta) - 3 t eps
    a2 = 3 ( t (1+2 alpha+beta) - beta )
    a1 = t (-(3+3 alpha+beta) + 3 u) + 3 (alpha+beta) + 3 t eps
    a0 = -(1+2 alpha+beta)

The eps-dependence is EXACTLY linear: P(y;u,eps) = P(y;u,0) + 3 t eps (y-y^3).
The instability boundary u*(eps) is the first positive-to-negative zero of the
universal X(u); the branch constants are

    c^2  = sin^2 theta [u* - G(delta*)],   G(delta) = delta(1+alpha delta+beta delta^2/3)
    kappa = 1 + d/deps [u*(eps) - G(delta*(eps))]_{eps=0}

Könnig (2015, arXiv:1503.07436, Sec. III A): on the finite branch the
Higuchi/ghost bound is equivalent to dy/dN >= 0; a negative dy/dN anywhere on
the branch signals the ghost regime. This module re-derives everything from
the reduced quartic alone and is cross-checked against the full Model
(src/core.py + src/scalarsec.py) in the deep self-similar regime
(m_FP >> 1, eps << 1). Units H0 = 1.
"""
import numpy as np
from scipy.optimize import brentq

from core import Model
import scalarsec

_ZSCAN = np.geomspace(0.01, 1e4, 400)


def zstar_any(m, Om):
    """First X(z) > 0 -> X < 0 crossing on a wide geometric grid."""
    Xg = np.array([scalarsec.X_of_z(m, z, Om) for z in _ZSCAN])
    for i in range(len(_ZSCAN) - 1):
        if Xg[i] > 0.0 and Xg[i + 1] < 0.0:
            return brentq(lambda z: scalarsec.X_of_z(m, z, Om),
                          _ZSCAN[i], _ZSCAN[i + 1], xtol=1e-12, rtol=1e-14)
    raise ValueError("no z* crossing found")

Y_LO = 1e-12
U_SCAN = np.geomspace(1e-3, 4000.0, 400)
U_DEEP = (4e3, 1e4, 1e5, 1e6, 1e8)


def reduced_coeffs(t, alpha, beta, u, eps):
    """Reduced-quartic coefficients (a4, a3, a2, a1, a0), highest degree first."""
    a4 = t * beta
    a3 = -3.0 * t * (alpha + beta) + (1.0 - alpha + beta) - 3.0 * t * eps
    a2 = 3.0 * (t * (1.0 + 2.0 * alpha + beta) - beta)
    a1 = t * (-(3.0 + 3.0 * alpha + beta) + 3.0 * u) + 3.0 * (alpha + beta) + 3.0 * t * eps
    a0 = -(1.0 + 2.0 * alpha + beta)
    return np.array([a4, a3, a2, a1, a0], dtype=float)


def P_poly(coeffs, y):
    return ((((coeffs[0] * y + coeffs[1]) * y + coeffs[2]) * y + coeffs[3]) * y
            + coeffs[4])


def dP_dy(coeffs, y):
    return (4.0 * coeffs[0] * y ** 3 + 3.0 * coeffs[1] * y ** 2
            + 2.0 * coeffs[2] * y + coeffs[3])


def y_of_u(t, alpha, beta, u, eps):
    """Smallest positive real root of the reduced quartic (finite branch)."""
    coeffs = reduced_coeffs(t, alpha, beta, u, eps)
    if coeffs[4] == 0.0:
        raise ValueError("degenerate branch (1+2*alpha+beta = 0): no finite branch")
    f = lambda y: P_poly(coeffs, y)
    if f(Y_LO) * f(1.0) > 0.0:
        raise ValueError("no sign change on (0, 1): no finite-branch root")
    return brentq(f, Y_LO, 1.0, xtol=1e-14, rtol=1e-14)


def dydN_of_u(t, alpha, beta, u, eps):
    """dy/dN on the reduced problem, N = ln a (u ~ a^{-3})."""
    y = y_of_u(t, alpha, beta, u, eps)
    coeffs = reduced_coeffs(t, alpha, beta, u, eps)
    # dP/du = 3 t y (the a1*y term), so dy/du = -3 t y/(dP/dy);
    # dy/dN = -3u dy/du = 9 t u y/(dP/dy)
    return 9.0 * t * u * y / dP_dy(coeffs, y)


def G_of(delta, alpha, beta):
    return delta * (1.0 + alpha * delta + beta / 3.0 * delta ** 2)


def X_of_u(t, alpha, beta, u, eps):
    """Universal X(u) on the reduced problem (Akrami et al. 2015, Eq. 38;
    B_n/(s^2 m_FP^2) with Omega_r = 0)."""
    y = y_of_u(t, alpha, beta, u, eps)
    yp = dydN_of_u(t, alpha, beta, u, eps)
    b1 = 1.0 + 2.0 * alpha + beta
    b2 = -(alpha + beta)
    b3 = beta
    T1 = ((b1 + 4.0 * b2 * y + 3.0 * b3 * y ** 2) * yp
          / (3.0 * y * (b1 + 2.0 * b2 * y + b3 * y ** 2)))
    T2 = ((1.0 + t * y ** 2) * (b1 - b3 * y ** 2) * yp ** 2
          / (3.0 * t * y ** 3 * 3.0 * u))
    return 1.0 + T1 - T2


def find_ustar(t, alpha, beta, eps, u_hi=4000.0):
    """First positive-to-negative zero of X(u) on the reduced problem."""
    Xu = np.array([X_of_u(t, alpha, beta, u, eps) for u in U_SCAN])
    for i in range(len(U_SCAN) - 1):
        if Xu[i] > 0.0 and Xu[i + 1] < 0.0:
            return brentq(lambda u: X_of_u(t, alpha, beta, u, eps),
                          U_SCAN[i], U_SCAN[i + 1], xtol=1e-12, rtol=1e-14)
    raise ValueError("no X(u) = 0 crossing found below u = %g" % u_hi)


def F_of_eps(t, alpha, beta, eps):
    """The combination entering H^2(z*) scaled by 1/(s^2 m_FP^2):
    F(eps) = u*(eps) - G(delta*(eps))."""
    u = find_ustar(t, alpha, beta, eps)
    y = y_of_u(t, alpha, beta, u, eps)
    return u - G_of(1.0 - y, alpha, beta)


def branch_constants(theta_deg=6.1, alpha=0.0, beta=3.0, h_eps=1e-5):
    """Branch constants (u*, y*, c, kappa) from the reduced problem alone."""
    th = np.radians(theta_deg)
    t = np.tan(th) ** 2
    s2 = np.sin(th) ** 2
    u0 = find_ustar(t, alpha, beta, 0.0)
    y0 = y_of_u(t, alpha, beta, u0, 0.0)
    c = np.sqrt(s2 * (u0 - G_of(1.0 - y0, alpha, beta)))
    kappa = 1.0 + (F_of_eps(t, alpha, beta, h_eps)
                   - F_of_eps(t, alpha, beta, -h_eps)) / (2.0 * h_eps)
    return dict(theta_deg=theta_deg, alpha=alpha, beta=beta,
                tan2=t, sin2=s2, ustar=u0, ystar=y0, c=c, kappa=kappa)


def u_of_y(theta_deg, y, alpha=0.0, beta=3.0):
    """Explicit rational inverse of the reduced quartic at eps = 0, benchmark
    branch (alpha=0, beta=3): solving P(y; u, 0) = 0 for u (u enters only
    through a1)."""
    th = np.radians(theta_deg)
    t = np.tan(th) ** 2
    return ((1.0 - y) * (3.0 * t * y ** 3 - 6.0 * t * y ** 2 + 6.0 * t * y
                         + 4.0 * y ** 2 - 5.0 * y + 4.0) / (3.0 * t * y))


def q_coeffs(theta_deg, alpha=0.0, beta=3.0):
    """DESCENDING-power coefficients (np.polyval/np.roots convention) of the
    degree-10 numerator polynomial Q(y) of the exact rational boundary
    function X(y) = -Q(y)/[(3y^2-6y+4)A(y)^2] (benchmark branch alpha=0,
    beta=3, eps=0); A(y) = 9ty^4-18ty^3+12ty^2+8y^3-9y^2+4. Both denominator
    factors are strictly positive on 0<y<1 (gate-checked), so the zeros of X
    on the branch are exactly the roots of Q. Caution: the natural reading
    order of the formula is ascending in y — do not pass this list to
    polyval/roots without keeping the descending convention."""
    th = np.radians(theta_deg)
    t = np.tan(th) ** 2
    return [81.0 * t ** 2,
            216.0 * t - 486.0 * t ** 2,
            96.0 - 1152.0 * t + 1404.0 * t ** 2,
            -432.0 + 2976.0 * t - 2376.0 * t ** 2,
            979.0 - 4842.0 * t + 2376.0 * t ** 2,
            -1542.0 + 5340.0 * t - 1296.0 * t ** 2,
            1746.0 - 4140.0 * t + 288.0 * t ** 2,
            -1376.0 + 2304.0 * t,
            744.0 - 912.0 * t,
            -288.0 + 192.0 * t,
            64.0]

def X_rational(theta_deg, y, alpha=0.0, beta=3.0):
    t = np.tan(np.radians(theta_deg)) ** 2
    q = np.polyval(np.array(q_coeffs(theta_deg, alpha, beta)), y)
    A = 9.0 * t * y ** 4 - 18.0 * t * y ** 3 + 12.0 * t * y ** 2 + 8.0 * y ** 3 - 9.0 * y ** 2 + 4.0
    return -q / ((3.0 * y ** 2 - 6.0 * y + 4.0) * A ** 2)


def sturm_root_count(c, x_lo, x_hi):
    """Number of real roots of the polynomial (descending coefficients c)
    in (x_lo, x_hi) via the Sturm sequence (numeric evaluation)."""
    def divmod_poly(p, q):
        p = [float(v) for v in p]
        out = [0.0] * (len(p) - len(q) + 1)
        while len(p) >= len(q) and any(abs(v) > 0 for v in p):
            f = p[0] / q[0]
            out[len(out) - (len(p) - len(q)) - 1] = f
            p = [pc - f * qc for pc, qc in zip(p, q + [0.0] * (len(p) - len(q)))]
            p = p[1:]
        return out, p

    def var_count(seq, x):
        vals = []
        for poly in seq:
            v = 0.0
            for cf in poly:
                v = v * x + cf
            vals.append(v)
        vals = [v for v in vals if abs(v) > 1e-10]
        return sum(1 for i in range(len(vals) - 1) if vals[i] * vals[i + 1] < 0)

    seq = [[float(v) for v in c], np.polyder(np.array(c)).tolist()]
    while len(seq[-1]) > 1:
        _, r = divmod_poly(seq[-2], seq[-1])
        r = [-v for v in r]
        while len(r) > 1 and abs(r[0]) < 1e-14:
            r = r[1:]
        seq.append(r)
        if len(seq) > 20:
            break
    return var_count(seq, x_lo) - var_count(seq, x_hi)


def X_asymptote(theta_deg=6.1, alpha=0.0, beta=3.0, u_deep=1e8):
    """Analytic asymptote of the universal X(u) as u -> infinity.
    Dominant balance a1 y + a0 = 0 gives y ~ (1+2a+b)/(3 t u); then
    dy/dN = 3 y, T1 -> 1, T2 -> 3, hence X = 1 + T1 - T2 -> -1 exactly."""
    th = np.radians(theta_deg)
    t = np.tan(th) ** 2
    y = y_of_u(t, alpha, beta, u_deep, 0.0)
    X = X_of_u(t, alpha, beta, u_deep, 0.0)
    y_asym = (1.0 + 2.0 * alpha + beta) / (3.0 * t * u_deep)
    return dict(u=u_deep, y=y, X=X, y_asym=y_asym, X_asym=-1.0)


def X_deep_negative(theta_deg=6.1, alpha=0.0, beta=3.0):
    """X(u) < 0 at all deep-grid points (beyond the u* scan range)."""
    th = np.radians(theta_deg)
    t = np.tan(th) ** 2
    return [(u, X_of_u(t, alpha, beta, u, 0.0)) for u in U_DEEP]


def zstar_master(theta_deg=6.1, alpha=0.0, beta=3.0, Omega_L=0.699,
                 Omega_r=9.08e-5, mFP_H0=25.0):
    """Master-equation prediction of z*: Omega_m(1+z*)^3 =
    sin^2(theta) m_FP^2 u*(eps) [1 + O(rho)], with the self-consistent
    closure Omega_m and the exact reduced-problem u*(eps)."""
    th = np.radians(theta_deg)
    t = np.tan(th) ** 2
    s2 = np.sin(th) ** 2
    m = Model(theta_deg=theta_deg, mFP_H0=float(mFP_H0), alpha=alpha, beta=beta,
              Omega_L=Omega_L, Omega_r=Omega_r)
    Om = m.solve_Om0()
    eps = Omega_L / (s2 * m.mFP2)
    u = find_ustar(t, alpha, beta, eps)
    return (u * s2 * m.mFP2 / Om) ** (1.0 / 3.0) - 1.0


def zstar_master_check(theta_deg=6.1, alpha=0.0, beta=3.0, Omega_L=0.699,
                       Omega_r=9.08e-5,
                       m_grid=(3.0, 6.0, 10.0, 15.0, 25.0, 40.0, 60.0, 100.0, 200.0)):
    """Master equation vs the full Model across the self-similar regime.
    Residual is O(rho(z*)) -- purely radiative once u*(eps) is included."""
    th = np.radians(theta_deg)
    s2 = np.sin(th) ** 2
    out = []
    for mf in m_grid:
        m = Model(theta_deg=theta_deg, mFP_H0=float(mf), alpha=alpha, beta=beta,
                  Omega_L=Omega_L, Omega_r=Omega_r)
        Om = m.solve_Om0()
        zg = np.geomspace(1e-3, 3000.0, 500)
        Xg = np.array([scalarsec.X_of_z(m, z, Om) for z in zg])
        z_full = None
        for i in range(len(zg) - 1):
            if Xg[i] > 0.0 and Xg[i + 1] < 0.0:
                z_full = brentq(lambda z: scalarsec.X_of_z(m, z, Om),
                                zg[i], zg[i + 1], xtol=1e-12)
                break
        z_m = zstar_master(theta_deg=theta_deg, alpha=alpha, beta=beta,
                           Omega_L=Omega_L, Omega_r=Omega_r, mFP_H0=float(mf))
        out.append(dict(mFP_H0=float(mf), eps=Omega_L / (s2 * m.mFP2),
                        zstar_full=z_full, zstar_master=z_m,
                        dev=abs(z_m - z_full) / z_full))
    return out


def eps_linearity_slope(theta_deg=6.1, alpha=0.0, beta=3.0, h_eps=1e-5):
    """d u*/deps at eps = 0 (central difference on the reduced problem)."""
    th = np.radians(theta_deg)
    t = np.tan(th) ** 2
    return (find_ustar(t, alpha, beta, h_eps)
            - find_ustar(t, alpha, beta, -h_eps)) / (2.0 * h_eps)


def eps_shift_coeffs(t, eps, alpha=0.3, beta=2.0, u=5.0):
    """Exact eps-shifts of the reduced-quartic coefficients (Lemma 1):
    da3 = -3 t eps, da1 = +3 t eps, all others zero."""
    c0 = reduced_coeffs(t, alpha, beta, u, 0.0)
    c1 = reduced_coeffs(t, alpha, beta, u, eps)
    return c1 - c0


def ghost_free_scan(theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0,
                    Omega_L=0.699, Omega_r=9.08e-5, n_grid=2000):
    """Min dy/dN over the observable branch 0 <= z <= 1090 (Könnig 2015:
    the Higuchi/ghost bound is dy/dN >= 0 on the finite branch)."""
    m = Model(theta_deg=theta_deg, mFP_H0=mFP_H0, alpha=alpha, beta=beta,
              Omega_L=Omega_L, Omega_r=Omega_r)
    Om = m.solve_Om0()
    zg = np.concatenate(([0.0], np.geomspace(1e-3, 1090.0, n_grid - 1)))
    yp = np.array([scalarsec.dy_dN(m, z, Om) for z in zg])
    return dict(z=zg, dydN=yp, min=float(yp.min()), z_at_min=float(zg[yp.argmin()]))


def full_model_crosscheck(theta_deg=6.1, alpha=0.0, beta=3.0, Omega_L=0.699,
                          mFP_H0=2000.0, Omega_r=0.0):
    """Cross-check the reduced problem against the full Model in the deep
    self-similar regime (m_FP = 2000, Omega_r = 0): locate the boundary in
    u-space via the full X(z), and compare H(z*)/m_FP with c."""
    m = Model(theta_deg=theta_deg, mFP_H0=mFP_H0, alpha=alpha, beta=beta,
              Omega_L=Omega_L, Omega_r=Omega_r)
    Om = m.solve_Om0()
    s2m2 = m.sin2 * m.mFP2

    def z_of_u(u):
        return (u * s2m2 / Om) ** (1.0 / 3.0) - 1.0

    u_full = brentq(lambda u: scalarsec.X_of_z(m, z_of_u(u), Om), 1.0, 200.0,
                    xtol=1e-10, rtol=1e-12)
    z_full = z_of_u(u_full)
    r_full = np.sqrt(m.E2_bim(z_full, Om)) / mFP_H0
    th = np.radians(theta_deg)
    ss = branch_constants(theta_deg=theta_deg, alpha=alpha, beta=beta)
    return dict(ustar_full=u_full, ustar_reduced=ss["ustar"],
                r_full=r_full, c=ss["c"],
                eps=Omega_L / s2m2)


def zstar_scaling(theta_deg=6.1, alpha=0.0, beta=3.0, Omega_L=0.699,
                  m_grid=(10.0, 20.0, 40.0, 80.0, 160.0), Omega_r=9.08e-5):
    """Slope of ln(1+z*) vs ln m_FP (self-similar prediction: 2/3)."""
    pts = []
    for mf in m_grid:
        m = Model(theta_deg=theta_deg, mFP_H0=float(mf), alpha=alpha, beta=beta,
                  Omega_L=Omega_L, Omega_r=Omega_r)
        Om = m.solve_Om0()
        zs = zstar_any(m, Om)
        pts.append((np.log(mf), np.log1p(zs)))
    lm = np.array([p[0] for p in pts])
    lz = np.array([p[1] for p in pts])
    return float(np.polyfit(lm, lz, 1)[0])


def plateau_scan(theta_deg=6.1, alpha=0.0, beta=3.0, Omega_L=0.699,
                 m_grid=(500.0, 1000.0, 2000.0, 4000.0, 8000.0)):
    """r = H(z*)/m_FP at Omega_r = 0, deep self-similar regime (plateau = c)."""
    out = []
    for mf in m_grid:
        m = Model(theta_deg=theta_deg, mFP_H0=float(mf), alpha=alpha, beta=beta,
                  Omega_L=Omega_L, Omega_r=0.0)
        Om = m.solve_Om0()
        zs = zstar_any(m, Om)
        out.append((mf, np.sqrt(m.E2_bim(zs, Om)) / mf, Omega_L / (m.sin2 * m.mFP2)))
    return out


def Hstar_law_check(theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0,
                    Omega_L=0.699, Omega_r=9.08e-5):
    """H(z*) from the full model vs sqrt(c^2 m_FP^2 + kappa Omega_L)."""
    m = Model(theta_deg=theta_deg, mFP_H0=mFP_H0, alpha=alpha, beta=beta,
              Omega_L=Omega_L, Omega_r=Omega_r)
    Om = m.solve_Om0()
    zs = scalarsec.find_zstar(m, Om)
    H_full = np.sqrt(m.E2_bim(zs, Om))
    ss = branch_constants(theta_deg=theta_deg, alpha=alpha, beta=beta)
    H_law = np.sqrt(ss["c"] ** 2 * mFP_H0 ** 2 + ss["kappa"] * Omega_L)
    return dict(zstar=zs, H_full=H_full, H_law=H_law, rel=(H_law / H_full - 1.0))


if __name__ == "__main__":
    ss = branch_constants()
    print("branch constants (reduced problem):")
    for k in ("ustar", "ystar", "c", "kappa"):
        print("  %s = %.8f" % (k, ss[k]))
    print("  d u*/deps (eps=0) = %.6f" % eps_linearity_slope())
    gf = ghost_free_scan()
    print("ghost-free scan: min dy/dN = %.6e at z = %.4f (> 0 required)" %
          (gf["min"], gf["z_at_min"]))
    xc = full_model_crosscheck()
    print("full-model crosscheck (m = 2000, Omega_r = 0): "
          "u* full %.6f vs reduced %.6f; r %.8f vs c %.8f (eps = %.2e)"
          % (xc["ustar_full"], xc["ustar_reduced"], xc["r_full"], xc["c"],
             xc["eps"]))
    print("z* scaling slope: %.4f" % zstar_scaling())
    print("plateau r = H(z*)/m_FP (Omega_r = 0):",
          ["%.8f" % r for _, r, _ in plateau_scan()])
    hl = Hstar_law_check()
    print("H(z*) law: full %.6f vs law %.6f (rel %+.4f), z* = %.6f"
          % (hl["H_full"], hl["H_law"], hl["rel"], hl["zstar"]))
    xa = X_asymptote()
    print("X(u) asymptote at u=1e8: X = %.8f (analytic -1), y = %.4e vs %.4e"
          % (xa["X"], xa["y"], xa["y_asym"]))
    print("X(u) deep grid:", ["%.4f" % X for _, X in X_deep_negative()])
    print("rational X(y): Sturm roots in (0,1) =",
          sturm_root_count(q_coeffs(6.1), 1e-12, 1.0 - 1e-12),
          "; X_rational == X_of_u to 1e-12:",
          max(abs(X_rational(6.1, y) - X_of_u(np.tan(np.radians(6.1))**2, 0.0, 3.0, u_of_y(6.1, y), 0.0))
              for y in (0.3, 0.5, 0.9)))
    print("master equation vs full model:")
    for r in zstar_master_check():
        print("  m=%6.1f eps=%7.4f  z* full=%9.4f  master=%9.4f  dev=%.3f%%"
              % (r["mFP_H0"], r["eps"], r["zstar_full"], r["zstar_master"],
                 100.0 * r["dev"]))
