"""Paper figures (Figs. 1-7), regenerated from the code base.

Every curve is either computed here with the canonical modules (core.py,
scalarsec.py, selfsim.py) or read from results/*.json written by run_all.sh.
Numbers that appear on the figures or in their captions are asserted against
the values quoted in the paper (results/figure_checks.txt); the script exits
non-zero if any of them moves.

    python3 make_figures.py            # all figures
    python3 make_figures.py 1 4 7      # selected figures
"""
import sys, json, pathlib
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter, LogLocator

SRC = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SRC))
from core import Model, E2_lcdm
import scalarsec
import selfsim
import figstyle as fs
from figstyle import (INK, INK2, BLUE, ORANGE, AQUA, VIOLET, RED, GREY,
                      panel_label, save, RESDIR)

CK = fs.Checks()


def load(name):
    return json.load(open(RESDIR / name))


def benchmark():
    m = Model()                       # theta=6.1, m_FP/H0=25, (alpha,beta)=(0,3), Omega_L=0.699
    return m, m.solve_Om0()


def Omega_DE(m, z, Om):
    d = 1.0 - m.solve_y(z, Om)
    return m.Omega_L - m.sin2 * m.mFP2 * d * (1 + m.alpha * d + m.beta / 3.0 * d ** 2)


def z_screen(m, Om):
    """Vainshtein boundary H(z_V) = m_FP (order-of-magnitude criterion)."""
    return brentq(lambda z: np.sqrt(m.E2_bim(z, Om)) - m.mFP_H0, 0.05, 3000)


def sci(x, pos=None):
    return r"$10^{%d}$" % int(round(np.log10(x)))


# ---------------------------------------------------------------- Fig. 1
def fig1_background():
    m, Om = benchmark()
    z0 = brentq(lambda z: Omega_DE(m, z, Om), 3.0, 6.0, xtol=1e-12)
    zV = z_screen(m, Om)
    CK.close("Fig1 Omega_m closure", Om, 0.3043547, 1e-7)
    CK.close("Fig1 z0 (Omega_DE zero crossing)", z0, 4.882, 6e-4)
    CK.close("Fig1 z_V (H = m_FP)", zV, 11.72, 6e-3)

    zz = np.linspace(0.0, 12.0, 481)
    y = np.array([m.solve_y(z, Om) for z in zz])
    ode = np.array([Omega_DE(m, z, Om) for z in zz])
    e2b = np.array([m.E2_bim(z, Om) for z in zz])
    OmL = 1.0 - Om - m.Omega_r
    e2l = E2_lcdm(zz, Om, m.Omega_r)
    dE = 100.0 * (e2b / e2l - 1.0)

    # Table A2 / A3 anchor values
    def at(z): return int(np.argmin(abs(zz - z)))
    CK.close("Fig1 y(12)", y[at(12.0)], 0.442652, 2e-6)
    CK.close("Fig1 dE2(z=5) [%]", dE[at(5.0)], -1.1088, 2e-4)
    CK.close("Fig1 dE2(z=12) [%]", dE[at(12.0)], -0.7667, 2e-4)
    CK.close("Fig1 max |dE2| [%] (z<=12)", abs(dE).max(), 1.11, 6e-3)
    CK.close("Fig1 z at max |dE2|", zz[int(np.argmax(abs(dE)))], 5.0, 0.4)
    CK.close("Fig1 Omega_DE(inf) asymptote",
             m.Omega_L - m.sin2 * m.mFP2 * (1 + m.alpha + m.beta / 3.0), -13.42, 5e-3)
    CK.close("Fig1 w_DE(0)", m.w_DE(0.0, Om), -1.0050, 1e-4)
    CK.close("Fig1 w_DE(4)", m.w_DE(4.0, Om), -2.5988, 1e-3)

    zw = np.linspace(0.0, 4.75, 191)
    w = np.array([m.w_DE(z, Om) for z in zw])

    fig, ax = plt.subplots(2, 2, figsize=(6.5, 4.9), constrained_layout=True)
    a, b, c, d = ax.ravel()
    for axx in (a, b, c):
        axx.axvline(zV, color=GREY, lw=0.9, ls=(0, (4, 3)))
        axx.axvline(z0, color=INK2, lw=0.9, ls=(0, (1, 2)))
    # (a) y(z)
    a.plot(zz, y, color=BLUE)
    a.set_xlabel(r"redshift $z$"); a.set_ylabel(r"$y(z)=a_f/a_g$")
    a.set_ylim(0.38, 1.02)
    a.text(zV - 0.25, 0.97, r"$z_V$", color=INK2, ha="right", va="top")
    a.text(z0 - 0.25, 0.97, r"$z_0$", color=INK2, ha="right", va="top")
    a.text(0.4, 0.45, r"de Sitter attractor: $y\to1$ as $z\to-1$", color=INK2, ha="left")
    panel_label(a, "(a)", x=-0.02, y=1.10)
    # (b) Omega_DE(z)
    b.plot(zz, ode, color=ORANGE)
    b.axhline(0, color=INK2, lw=0.7)
    b.axhline(-13.42, color=INK2, lw=0.8, ls=(0, (6, 3)))
    b.text(0.4, -12.6, r"asymptote $-13.42$", color=INK2, va="bottom")
    b.plot([z0], [0], "o", color=ORANGE, mec="white", mew=0.8, ms=6, zorder=5)
    b.annotate(r"$z_0=4.882$", xy=(z0, 0), xytext=(5.6, 3.2), color=INK,
               arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
    b.set_ylim(-14.5, 6.0)
    b.set_xlabel(r"redshift $z$"); b.set_ylabel(r"$\Omega_{\rm DE}(z)$")
    panel_label(b, "(b)", x=-0.02, y=1.10)
    # (c) dE2
    c.plot(zz, dE, color=VIOLET)
    c.axhline(0, color=INK2, lw=0.7)
    ib = int(np.argmax(abs(dE)))
    c.plot([zz[ib]], [dE[ib]], "o", color=VIOLET, mec="white", mew=0.8, ms=6, zorder=5)
    c.annotate(r"max $|\delta E^2|=1.11\%$" + "\n" + r"near $z=5$", xy=(zz[ib], dE[ib]),
               xytext=(5.3, -0.28), color=INK, ha="left",
               arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
    c.set_xlabel(r"redshift $z$")
    c.set_ylabel(r"$\delta E^2$ [%]")
    c.set_ylim(-1.35, 0.1)
    panel_label(c, "(c)", x=-0.02, y=1.10)
    # (d) w_DE
    d.axhspan(-6, -1, color="#f1d9cf", alpha=0.55, lw=0)
    d.axhline(-1, color=INK2, lw=0.8)
    d.plot(zw, w, color=AQUA)
    d.set_ylim(-4.2, -0.8); d.set_xlim(0, 4.9)
    d.axvline(z0, color=GREY, lw=0.9, ls=(0, (1, 2)))
    d.text(z0 - 0.08, -1.12, r"$w_{\rm DE}\to-\infty$ at $z_0$", ha="right", va="top", color=INK2)
    d.text(0.1, -3.9, "phantom, $w_{\\rm DE}<-1$", color=INK2, va="bottom")
    d.set_xlabel(r"redshift $z$"); d.set_ylabel(r"$w_{\rm DE}(z)$")
    panel_label(d, "(d)", x=-0.02, y=1.10)
    save(fig, "fig1_background")


# ---------------------------------------------------------------- Fig. 2
def fig2_mechanism():
    m, Om = benchmark()
    Om_eff = m.cos2 * Om
    CK.close("Fig2 Omega_eff = cos^2(theta) Omega_m (BAO closure 0.30439)",
             m.cos2 * 0.30439172, 0.30096, 3e-5)
    zz = np.concatenate(([0.0], np.geomspace(1e-2, 12.0, 400)))
    dev = np.array([m.E2_bim(z, Om) / (m.Omega_L + Om_eff * (1 + z) ** 3
                                       + m.Omega_r * (1 + z) ** 4) - 1.0 for z in zz])
    win = zz <= 2.33
    CK.true("Fig2 dark degeneracy: max |dev| for z<=2.33 below 4e-5 (paper 3.4e-5, gate): %.2e"
            % abs(dev[win]).max(), abs(dev[win]).max() < 4e-5)

    i5, i12 = int(np.argmin(abs(zz - 5.0))), int(np.argmin(abs(zz - 12.0)))
    CK.close("Fig2 residual at z=5 (text: 1.6e-5)", dev[i5], 1.58e-5, 1e-7)
    CK.close("Fig2 residual at z=12 (paper 3.6e-3)", dev[i12], 3.6e-3, 6e-5)
    CK.true("Fig2 residual changes sign once, near z = 5", int(np.sum(np.diff(np.sign(dev[1:])) != 0)) == 1)

    S = load("summary.json")
    bao, sd, sf = S["bao"], S["sne"]["diag"], S["sne"]["full"]
    CK.close("Fig2 BAO dchi2", bao["dchi2"], 0.179, 5e-4)
    CK.close("Fig2 BAO density-shift term", bao["Om_shift_penalty"], 0.665, 5e-4)
    CK.close("Fig2 BAO E^2-difference term", bao["bim_correction"], -0.486, 5e-4)
    CK.close("Fig2 SNe full dchi2", sf["dchi2"], 2.94, 5e-3)
    CK.close("Fig2 SNe full pure", sf["pure"], 0.64, 5e-3)
    CK.close("Fig2 SNe diag dchi2", sd["dchi2"], 7.32, 5e-3)
    CK.close("Fig2 SNe diag pure", sd["pure"], 1.01, 5e-3)
    cmb0 = S["cmb_scan"]["6.1"]["dchi2"]
    CK.close("Fig2 CMB dchi2", cmb0, 0.0001702627819509077, 5e-5)

    fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 2.75), constrained_layout=True,
                               gridspec_kw=dict(width_ratios=[1.12, 1]))
    # (a) residual to the effective law
    a.axvspan(0, 2.33, color="#dbe7f7", alpha=0.8, lw=0)
    a.text(1.05, 0.9, "BAO + SNe\n" + r"$z\leq2.33$", transform=a.get_xaxis_transform(),
           ha="center", va="top", color=BLUE, fontsize=7.5)
    a.plot(zz, dev, color=VIOLET)
    a.set_yscale("symlog", linthresh=1e-6, linscale=0.4)
    a.set_xlim(0, 12); a.set_ylim(-1e-4, 1e-2)
    a.axhline(0, color=INK2, lw=0.7)
    a.axhline(3.4e-5, color=GREY, lw=0.8, ls=(0, (4, 3)))
    a.axhline(-3.4e-5, color=GREY, lw=0.8, ls=(0, (4, 3)))
    a.text(11.8, 4.4e-5, r"$3.4\times10^{-5}$", ha="right", va="bottom", color=INK2)
    a.text(11.8, -3.4e-5, r"$-3.4\times10^{-5}$", ha="right", va="top", color=INK2)
    a.set_yticks([-1e-4, -1e-5, 0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
    a.set_yticklabels([r"$-10^{-4}$", r"$-10^{-5}$", "0", r"$10^{-6}$", r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$"])
    a.set_xlabel(r"redshift $z$")
    a.set_ylabel(r"$E^2_{\rm bim}/E^2_{\rm eff}-1$")
    panel_label(a, "(a)", x=-0.02, y=1.10)
    # (b) decomposition of the per-probe delta chi^2
    labels = ["BAO", "SNe\n(full cov.)", "SNe\n(diag.)", "CMB"]
    shift = [bao["Om_shift_penalty"], sf["dchi2"] - sf["pure"], sd["dchi2"] - sd["pure"], 0.0]
    corr = [bao["bim_correction"], sf["pure"], sd["pure"], cmb0]
    tot = [bao["dchi2"], sf["dchi2"], sd["dchi2"], cmb0]
    xs = np.arange(4); wbar = 0.6
    b.axhline(0, color=INK2, lw=0.7)
    b.bar(xs, shift, wbar, color=ORANGE, alpha=0.92, lw=0, label=r"density mismatch ($\Omega_{\rm m}$)")
    b.bar(xs, corr, wbar, bottom=[s if c > 0 else 0 for s, c in zip(shift, corr)],
          color=BLUE, alpha=0.92, lw=0, label=r"$E^2(z)$ correction")
    for x, t in zip(xs, tot):
        b.plot([x - 0.38, x + 0.38], [t, t], color=INK, lw=1.4)
        b.text(x + 0.42, t, r"$%s$" % ("%.3g" % t if t >= 0.01 else "%.4f" % t),
               ha="left", va="center", fontsize=7.5)
    b.set_xticks(xs); b.set_xticklabels(labels)
    b.set_xlim(-0.55, 3.95)
    b.set_ylabel(r"$\Delta\chi^2$ against $\Lambda$CDM")
    b.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0))
    b.set_ylim(-1.0, 9.6)
    b.grid(axis="x", visible=False)
    b.text(0.02, 0.03, "black bar: total", transform=b.transAxes, fontsize=7, color=INK2)
    panel_label(b, "(b)", x=-0.02, y=1.10)
    save(fig, "fig2_mechanism")



# ---------------------------------------------------------------- Fig. 3
def _quad_vertex(x, y):
    a, b, c = np.polyfit(x, y, 2)
    xv = -b / (2 * a)
    return xv, c - b * b / (4 * a)


def _cross(x, y, level, lo, hi):
    """Interpolated crossing of y(x) = level on [lo, hi] (cubic through the fine profile)."""
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(x, y)
    return brentq(lambda t: cs(t) - level, lo, hi)


def fig3_anchor_profile():
    A = load("anchor_scan.json")
    rows = sorted(A["profile_1D"] + A["profile_1D_fine"], key=lambda r: r[0])
    rows = np.array(rows)
    fine = np.array(sorted(A["profile_1D_fine"], key=lambda r: r[0]))
    coarse = np.array(A["profile_1D"])
    OL = rows[:, 0]
    # Table 6 checks
    def row(ol): return coarse[np.argmin(abs(coarse[:, 0] - ol))]
    for ol, want in ((0.600, (100.57, 12.95, 150.50, 264.02)), (0.699, (0.18, 2.94, 0.00, 3.12)),
                     (0.740, (21.78, 17.14, 42.04, 80.96)), (0.780, (111.71, 44.43, 195.29, 351.43))):
        r = row(ol)
        for nm, g, w in zip(("BAO", "SNe", "CMB", "joint"), r[1:], want):
            CK.close("Fig3 Table 6 %s at Omega_L=%.3f" % (nm, ol), g, w, 6e-3)
    # SNe-preferred anchor and the competing pulls
    sel = (fine[:, 0] >= 0.656) & (fine[:, 0] <= 0.684)
    xs_sne, ys_sne = _quad_vertex(fine[sel, 0], fine[sel, 2])
    CK.close("Fig3 SNe-preferred Omega_Lambda (paper ~0.668)", xs_sne, 0.668, 3e-3)
    r668 = np.array([np.interp(0.668, fine[:, 0], fine[:, k]) for k in (1, 2, 3, 4)])
    CK.close("Fig3 BAO penalty at 0.668 (paper +14.28)", r668[0], 14.28, 0.3)
    CK.close("Fig3 CMB penalty at 0.668 (paper +18.19)", r668[2], 18.19, 0.01)
    CK.close("Fig3 joint at 0.668 (paper 32.48)", r668[3], 32.48, 0.01)
    # joint vertex and tolerances (fine profile, 0.004 spacing)
    w = (fine[:, 0] >= 0.684) & (fine[:, 0] <= 0.712)
    xv, yv = _quad_vertex(fine[w, 0], fine[w, 4])
    CK.close("Fig3 joint vertex Omega_Lambda* (paper 0.6976)", xv, 0.6976, 1.5e-3)
    CK.close("Fig3 joint at vertex (paper 3.046)", yv, 3.046, 0.03)
    lo68 = _cross(fine[:, 0], fine[:, 4], yv + 1.0, 0.66, xv)
    hi68 = _cross(fine[:, 0], fine[:, 4], yv + 1.0, xv, 0.74)
    lo95 = _cross(fine[:, 0], fine[:, 4], yv + 3.84, 0.66, xv)
    hi95 = _cross(fine[:, 0], fine[:, 4], yv + 3.84, xv, 0.74)
    CK.close("Fig3 68%% lower (paper 0.692)", lo68, 0.692, 1.5e-3)
    CK.close("Fig3 68%% upper (paper 0.703)", hi68, 0.703, 1.5e-3)
    CK.close("Fig3 95%% lower (paper 0.687)", lo95, 0.687, 1.5e-3)
    CK.close("Fig3 95%% upper (paper 0.708)", hi95, 0.708, 1.5e-3)

    fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 2.9), constrained_layout=True,
                               gridspec_kw=dict(width_ratios=[1.25, 1]))
    series = [("BAO", 1, BLUE, "o"), ("SNe", 2, ORANGE, "s"), ("CMB", 3, AQUA, "^"),
              ("joint", 4, VIOLET, "D")]
    for nm, k, col, mk in series:
        a.plot(OL, rows[:, k], color=col, lw=1.5 if nm != "joint" else 2.0,
               ls="-" if nm != "joint" else (0, (5, 1.5)))
        a.plot(coarse[:, 0], coarse[:, k], mk, color=col, ms=3.6, mec="white", mew=0.5)
    a.set_yscale("symlog", linthresh=1.0, linscale=0.55)
    a.set_ylim(0, 450)
    a.set_yticks([0, 1, 10, 100]); a.set_yticklabels(["0", "1", "10", "100"])
    a.axvline(0.699, color=GREY, lw=0.9, ls=(0, (1, 2)))
    a.text(0.7015, 330, "benchmark\n" + r"$\Omega_\Lambda=0.699$", fontsize=7, color=INK2, va="top")
    a.plot([xs_sne], [0.0], "v", color=ORANGE, ms=5.5, mec="white", mew=0.5, clip_on=False, zorder=6)
    a.text(0.6, 3.2, "BAO", color=BLUE, fontsize=8, ha="left"); a.text(0.622, 12, "CMB", color=AQUA, ha="left")
    a.text(0.745, 9.0, "SNe", color=ORANGE, ha="left", va="center")
    a.text(0.705, 55, "joint", color=VIOLET, ha="left")
    a.set_xlim(0.595, 0.785)
    a.set_xlabel(r"anchor $\Omega_\Lambda$ at $\theta=6.1^\circ$")
    a.set_ylabel(r"$\Delta\chi^2$ against $\Lambda$CDM")
    panel_label(a, "(a)", x=-0.02, y=1.10)
    # (b) zoom on the joint profile, tolerance bands
    b.axvspan(lo95, hi95, color="#e4def7", alpha=0.9, lw=0)
    b.axvspan(lo68, hi68, color="#c9bff0", alpha=0.9, lw=0)
    b.plot(fine[:, 0], fine[:, 4], "-", color=VIOLET, lw=1.8)
    b.plot(fine[:, 0], fine[:, 4], "D", color=VIOLET, ms=3.2, mec="white", mew=0.4)
    b.plot([xv], [yv], "o", mfc="white", mec=INK, ms=5, zorder=6)
    b.plot([0.699], [3.1152], "*", color=INK, ms=8, zorder=7)
    b.annotate(r"vertex $0.6976$" + "\n" + r"(joint $3.05$)", xy=(xv, yv), xytext=(0.673, 9.5),
               fontsize=7.5, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7), ha="center")
    b.text(0.7035, 13.3, "$\\Delta\\chi^2$=1", fontsize=7, color="#4a3fa8", ha="left")
    b.text(0.7085, 11.3, "$\\Delta\\chi^2$=3.84", fontsize=7, color="#4a3fa8", ha="left")
    b.set_xlim(0.66, 0.74); b.set_ylim(0, 15)
    b.set_xlabel(r"anchor $\Omega_\Lambda$")
    b.set_ylabel(r"joint $\Delta\chi^2$")
    panel_label(b, "(b)", x=-0.02, y=1.10)
    save(fig, "fig3_anchor_profile")


# ---------------------------------------------------------------- Fig. 4
def fig4_anchor_plane():
    from scipy.interpolate import PchipInterpolator
    A = load("anchor_scan.json")
    grid = A["grid_2D"] + A["grid_2D_ext"]
    valley = A["theta_profile_continuous"]["rows"]       # [theta, OL*, joint]
    th_list = [4.0, 6.1, 8.0, 12.2, 18.3, 25.0, 30.0]
    V = {r[0]: (r[1], r[2]) for r in valley}
    # Table 7 / Sec. 7.4.3 checks
    for th, (ol, jj) in ((4.0, (0.6933, 3.15)), (5.0, (0.6951, 2.98)), (6.1, (0.6976, 3.04)),
                         (8.0, (0.7030, 4.41)), (12.2, (0.7197, 19.89)), (18.3, (0.7540, 114.64)),
                         (25.0, (0.8014, 435.72)), (30.0, (0.8378, 945.14))):
        CK.close("Fig4 Table 7 Omega_L* at %.1f deg" % th, V[th][0], ol, 6e-5)
        CK.close("Fig4 Table 7 joint min at %.1f deg" % th, V[th][1], jj, 6e-3)
    jmin = min(r[2] for r in valley)
    CK.close("Fig4 profile minimum (near 5 deg)", jmin, 2.98, 6e-3)
    CK.close("Fig4 GR limit (0.5 deg) above minimum", V[0.5][1] - jmin, 0.9, 0.06)

    fig = plt.figure(figsize=(6.5, 4.3), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.35, 1])
    a = fig.add_subplot(gs[:, 0]); b = fig.add_subplot(gs[0, 1]); c = fig.add_subplot(gs[1, 1])
    cm = plt.get_cmap("viridis")
    for i, th in enumerate(th_list):
        pts = {}
        for g in grid:
            if abs(g[0] - th) < 1e-9 and g[-1] is not None:
                pts[g[1]] = g[-1]
        pts[V[th][0]] = V[th][1]
        xs = np.array(sorted(pts)); ys = np.array([pts[x] for x in xs])
        xx = np.linspace(xs[0], xs[-1], 300)
        col = cm(0.08 + 0.80 * i / (len(th_list) - 1))
        a.plot(xx, 10 ** PchipInterpolator(xs, np.log10(ys))(xx), color=col, lw=1.4)
        a.plot(xs, ys, "o", color=col, ms=2.8, mec="white", mew=0.3)
        a.plot([V[th][0]], [V[th][1]], "o", mfc=col, mec=INK, mew=0.9, ms=5.5, zorder=6)
    ov = np.array([[r[1], r[2]] for r in valley])
    a.plot(ov[:, 0], ov[:, 1], color=INK, lw=0.8, ls=(0, (3, 2)), zorder=5)
    a.set_yscale("log"); a.set_ylim(2, 1800); a.set_xlim(0.59, 0.955)
    for th, dx, dy, ha in ((12.2, -0.008, 1.45, "right"), (18.3, -0.008, 1.5, "right"),
                           (25.0, -0.006, 1.22, "right"), (30.0, -0.006, 1.22, "right")):
        a.text(V[th][0] + dx, V[th][1] * dy, r"$%.1f^\circ$" % th, fontsize=7.5, ha=ha, va="bottom",
               color=cm(0.08 + 0.80 * th_list.index(th) / (len(th_list) - 1)))
    a.annotate(r"$4^\circ,\,6.1^\circ,\,8^\circ$", xy=(0.692, 3.25), xytext=(0.605, 2.6), fontsize=7.5,
               color=cm(0.2), ha="left", va="center",
               arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    a.set_xlabel(r"anchor $\Omega_\Lambda$")
    a.set_ylabel(r"joint $\Delta\chi^2$ (sum of per-probe values)")
    a.text(0.945, 2.3, "dots: grid points\ncircles: minimum at each $\\theta$\ndashed: valley (dark degeneracy)",
           fontsize=7, color=INK2, va="bottom", ha="right")
    panel_label(a, "(a)", x=-0.02, y=1.04)
    # (b) profile minimum vs theta
    th = np.array([r[0] for r in valley]); jm = np.array([r[2] for r in valley])
    b.axhspan(jmin, jmin + 1.0, color="#c9bff0", alpha=0.9, lw=0)
    b.plot(th, jm, "-o", color=VIOLET, ms=3.6, mec="white", mew=0.4)
    b.set_yscale("log"); b.set_ylim(2, 1300)
    b.set_xlabel(r"mixing angle $\theta$ [deg]")
    b.set_ylabel(r"min$_{\Omega_\Lambda}$ joint $\Delta\chi^2$")
    b.set_yticks([3, 10, 30, 100, 300, 1000]); b.set_yticklabels(["3", "10", "30", "100", "300", "1000"])
    b.yaxis.set_minor_locator(plt.NullLocator())
    panel_label(b, "(b)", x=-0.04, y=1.07)
    # (c) zoom
    z = th <= 8.0
    c.axhspan(jmin, jmin + 1.0, color="#c9bff0", alpha=0.9, lw=0)
    c.text(0.3, jmin + 0.5, r"$\Delta\chi^2\leq\min+1$", fontsize=7, va="center", color="#4a3fa8")
    c.plot(th[z], jm[z], "-o", color=VIOLET, ms=3.6, mec="white", mew=0.4)
    i61 = list(th).index(6.1)
    c.plot([6.1], [jm[i61]], "*", color=INK, ms=9, zorder=7)
    c.annotate("benchmark angle\n(profiled anchor)", xy=(6.1, jm[i61]), xytext=(5.3, 4.45), fontsize=7.5, ha="center",
               arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
    c.annotate("GR limit\n" + r"($\theta\to0$)", xy=(0.5, jm[0]), xytext=(1.9, 4.3), fontsize=7.5,
               arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
    c.set_xlim(0, 8.5); c.set_ylim(2.85, 4.7)
    c.set_xlabel(r"mixing angle $\theta$ [deg]")
    c.set_ylabel(r"joint $\Delta\chi^2$ (zoom)")
    panel_label(c, "(c)", x=-0.04, y=1.07)
    save(fig, "fig4_anchor_plane")


# ---------------------------------------------------------------- Fig. 5
def _lam_vertex(lams, chi2):
    lams = np.asarray(lams); chi2 = np.asarray(chi2)
    i = int(np.argmin(chi2)); i = min(max(i, 1), len(lams) - 2)
    a, b, c = np.polyfit(lams[i - 1:i + 2], chi2[i - 1:i + 2], 2)
    xv = -b / (2 * a)
    return xv, c - b * b / (4 * a)


def fig5_cmb_transplant():
    Sw = load("cmb_lambda_theta_sweep.json")
    ths = sorted(float(k) for k in Sw)
    cm = plt.get_cmap("viridis")
    # Table A5 checks
    for th, want in ((10.0, 0.48), (12.2, 0.58), (15.0, 0.65), (18.3, 0.69), (20.0, 0.71),
                     (25.0, 0.74), (30.0, 0.75)):
        d = Sw[str(th)]
        CK.close("Fig5 Table A5 lambda* at %.1f deg" % th, _lam_vertex(d["lambdas"], d["chi2"])[0],
                 want, 0.0075)
    CK.true("Fig5 lambda*(6.1 deg) < 0 (boundary)", _lam_vertex(Sw["6.1"]["lambdas"], Sw["6.1"]["chi2"])[0] < 0)
    CK.close("Fig5 chi2(6.1 deg, lambda=0)", Sw["6.1"]["chi2"][0], 0.0001702627819509077, 5e-5)
    CK.close("Fig5 chi2(6.1 deg, lambda=1)", Sw["6.1"]["chi2"][-1], 6.28707127758636, 6e-3)
    CK.close("Fig5 chi2(30 deg, lambda=0)", Sw["30.0"]["chi2"][0], 2412.984814454524, 6e-3)
    CK.close("Fig5 chi2(30 deg, lambda=0.75)", Sw["30.0"]["chi2"][3], 0.5554708132041974, 6e-3)
    CK.close("Fig5 chi2(30 deg, lambda=1)", Sw["30.0"]["chi2"][-1], 275.61452653073746, 6e-3)
    CK.close("Fig5 chi2(10 deg, lambda=0.5)", Sw["10.0"]["chi2"][2], 0.008640802322230948, 6e-4)
    best = {th: min(Sw[str(th)]["chi2"]) for th in ths}
    CK.true("Fig5 saturation: best chi2 on the lambda grid below ~2.5 (tolerance 2.55) for every theta >= 10 deg "
            "(max %.2f)" % max(v for t, v in best.items() if t >= 10.0),
            all(v < 2.55 for t, v in best.items() if t >= 10.0))

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(6.5, 2.6), constrained_layout=True,
                                  gridspec_kw=dict(width_ratios=[1.2, 1, 1.05]))
    sel = [6.1, 10.0, 15.0, 20.0, 30.0]
    for th in sel:
        d = Sw[str(th)]
        col = cm(0.08 + 0.8 * ths.index(th) / (len(ths) - 1))
        a.semilogy(d["lambdas"], np.maximum(d["chi2"], 2e-4), "-o", color=col, ms=3.4, mec="white", mew=0.4)
        xv, yv = _lam_vertex(d["lambdas"], d["chi2"])
        if 0 <= xv <= 1:
            a.plot([xv], [max(yv, 2e-4)], "v", color=col, ms=0, zorder=0)
        a.text(-0.04, max(d["chi2"][0], 2e-4), r"$%.1f^\circ$" % th, color=col, fontsize=7.5,
               va="center", ha="right")
    a.set_xlim(-0.30, 1.03); a.set_ylim(2e-4, 6e3)
    a.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    a.set_xlabel(r"transplant parameter $\lambda$")
    a.set_ylabel(r"compressed-CMB $\chi^2$")
    a.yaxis.set_major_locator(LogLocator(base=100.0, numticks=6))
    panel_label(a, "(a)", x=-0.02, y=1.10)
    # (b) lambda*(theta)
    tt = list(ths)
    lv = [_lam_vertex(Sw[str(t)]["lambdas"], Sw[str(t)]["chi2"])[0] for t in tt]
    b.axhspan(0, 1, color="#dbe7f7", alpha=0.8, lw=0)
    b.axhline(0, color=INK2, lw=0.7)
    b.text(12.5, 0.08, "admissible\n" + r"$\lambda^*\in[0,1]$", fontsize=7, color=BLUE, va="bottom")
    b.plot(tt[1:], lv[1:], "-o", color=VIOLET, ms=3.8, mec="white", mew=0.4)
    b.plot(tt[:2], lv[:2], "--", color=VIOLET, lw=1.0)
    b.plot(tt[:1], lv[:1], "o", mfc="white", mec=VIOLET, ms=4.6)
    b.set_xlim(4, 31.5); b.set_ylim(-0.32, 1.0)
    b.set_xlabel(r"mixing angle $\theta$ [deg]")
    b.set_ylabel(r"$\lambda^*(\theta)$")
    b.annotate(r"$6.1^\circ$: $\lambda^*<0$" + "\n(boundary)", xy=(tt[0], lv[0]), xytext=(8.4, -0.28),
               fontsize=7, color=INK2, va="bottom", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    panel_label(b, "(b)", x=-0.04, y=1.10)
    # (c) chi2 vs theta for fixed maps
    c0 = [Sw[str(t)]["chi2"][0] for t in ths]
    c1 = [Sw[str(t)]["chi2"][-1] for t in ths]
    cb = [best[t] for t in ths]
    c.semilogy(ths, np.maximum(c0, 2e-4), "-o", color=ORANGE, ms=3.4, mec="white", mew=0.4)
    c.semilogy(ths, c1, "-s", color=BLUE, ms=3.4, mec="white", mew=0.4)
    c.semilogy(ths, np.maximum(cb, 2e-4), "-^", color=AQUA, ms=3.8, mec="white", mew=0.4)
    c.set_ylim(2e-4, 6e3); c.set_xlim(5, 31)
    c.yaxis.set_major_locator(LogLocator(base=100.0, numticks=6))
    c.text(7.2, 1500, r"$\lambda=0$ (true density)", color=ORANGE, ha="left", va="center", fontsize=7.5)
    c.text(14.8, 4.2, r"$\lambda=1$ (effective)", color=BLUE, ha="left", va="bottom", fontsize=7.5, rotation=0)
    c.text(10.5, 0.0013, r"best of 5 $\lambda$ values", color=AQUA, ha="left", va="center", fontsize=7.5)
    c.set_xlabel(r"mixing angle $\theta$ [deg]")
    c.set_ylabel(r"compressed-CMB $\chi^2$")
    panel_label(c, "(c)", x=-0.04, y=1.10)
    save(fig, "fig5_cmb_transplant")



# ---------------------------------------------------------------- Fig. 6
def fig6_scalar_sector():
    m, Om = benchmark()
    zs = scalarsec.find_zstar(m, Om)
    zV = z_screen(m, Om)
    z0 = brentq(lambda z: Omega_DE(m, z, Om), 3.0, 6.0, xtol=1e-12)
    I1, _ = scalarsec.growth_integral(m, Om, zs)
    Iu, _ = quad(lambda z: np.sqrt(-scalarsec.X_of_z(m, z, Om)) / np.sqrt(m.E2_bim(z, Om)),
                 zs, zV, epsabs=1e-12, epsrel=1e-12, limit=300)
    frac = Iu / I1
    CK.close("Fig6 z*", zs, 7.8771, 1e-5)
    CK.close("Fig6 H(z*)/H0", np.sqrt(m.E2_bim(zs, Om)), 14.5573, 1e-4)
    CK.close("Fig6 I1", I1, 0.944870550, 1e-9)
    CK.close("Fig6 unscreened share of I1 (12.3%)", 100 * frac, 12.3, 0.05)
    for z, want in ((0.0, 0.999996), (5.0, 0.839634), (10.0, -0.466223), (1090.0, -1.163701)):
        CK.close("Fig6 X(%g)" % z, scalarsec.X_of_z(m, z, Om), want, 2e-6)
    CK.true("Fig6 ordering z0 < z* < z_V (%.3f < %.4f < %.3f)" % (z0, zs, zV), z0 < zs < zV)

    zg = np.geomspace(0.01, 1090.0, 700)
    X = np.array([scalarsec.X_of_z(m, z, Om) for z in zg])
    # cumulative growth integral from z* upward (trapezoid in ln z, dense grid)
    zc = np.geomspace(zs * (1 + 1e-9), 1090.0, 4000)
    f = np.array([np.sqrt(max(-scalarsec.X_of_z(m, z, Om), 0.0)) / np.sqrt(m.E2_bim(z, Om)) for z in zc])
    cum = cumulative_trapezoid(f * zc, np.log(zc), initial=0.0)
    CK.close("Fig6 cumulative I1 (trapezoid) vs quad", cum[-1], I1, 2e-3, rel=True)

    km = load("scalar_k_map.json")
    k_rec, k_g1, zh_g1 = km["k_at_zrec"], km["k_gamma1"], km["z_h_gamma1"]
    CK.close("Fig6 k(z_rec) = aH(z_rec) [Mpc^-1]", k_rec, 4.90e-3, 6e-3, rel=True)
    CK.close("Fig6 k for Gamma = 10 [Mpc^-1]", k_g1, 5.7e-4, 6e-3, rel=True)
    CK.close("Fig6 z_h(Gamma=10)", zh_g1, 18.5, 6e-3, rel=True)
    for r, (kk, zh, lg) in zip(km["rows"], ((0.001, 58.2, 1.76), (0.05, 20858, 87.87), (0.1, 43300, 175.74),
                                            (0.5, 2.23e5, 878.72), (1.0, 4.48e5, 1757.44), (5.0, 2.25e6, 8787.18))):
        CK.close("Fig6 Table 4 z_h(k=%g)" % kk, r["z_h"], zh, 6e-3, rel=True)
        CK.close("Fig6 Table 4 log10 Gamma(k=%g)" % kk, r["log10Gamma"], lg, 6e-3, rel=True)

    fig = plt.figure(figsize=(6.5, 5.0), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.05])
    a = fig.add_subplot(gs[0, :]); b = fig.add_subplot(gs[1, 0]); c = fig.add_subplot(gs[1, 1])
    # (a) X(z)
    a.axvspan(zs, 1090, color="#f6dccf", alpha=0.7, lw=0)
    a.axhline(0, color=INK2, lw=0.8)
    a.plot(zg, X, color=BLUE, lw=1.9)
    a.set_xscale("log"); a.set_xlim(0.01, 1090); a.set_ylim(-1.35, 1.2)
    for z, lab, dy in ((z0, r"$z_0=4.88$", 0.0), (zs, r"$z_*=7.877$", 0.0), (zV, r"$z_V=11.72$", 0.0)):
        a.axvline(z, color=GREY, lw=0.9, ls=(0, (3, 2)))
    a.text(z0 / 1.12, -1.05, r"$z_0$", ha="right", color=INK2); a.text(zs * 1.0, 1.0, r"$z_*$", ha="right", color=INK, va="center")
    a.text(zV * 1.12, -1.05, r"$z_V$", ha="left", color=INK2)
    a.text(30, 0.75, "gradient instability\n" + r"$X<0$ for $z>z_*=7.877$", color="#a2441f", ha="left", va="center")
    a.text(0.013, 0.55, r"$X(0)=+0.999996$", color=BLUE, ha="left")
    a.text(600, -1.31, r"$X(1090)=-1.1637$", color=BLUE, ha="right", fontsize=7.5)
    a.set_xlabel(r"redshift $z$"); a.set_ylabel(r"$X(z)=(aH/k)^2\omega^2$")
    panel_label(a, "(a)", x=-0.02, y=1.08)
    # (b) cumulative growth
    sc = zc >= zs
    b.axvspan(zs, zV, color="#c9e6d9", alpha=0.9, lw=0)
    b.plot(zc, cum / I1, color=ORANGE, lw=1.9)
    b.set_xscale("log"); b.set_xlim(zs, 1090); b.set_ylim(0, 1.02)
    b.axvline(zV, color=GREY, lw=0.9, ls=(0, (3, 2)))
    b.axhline(frac, color=GREY, lw=0.8, ls=(0, (1, 2)))
    b.text(zV * 1.1, frac - 0.03, r"$12.3\%$ at $z_V$", ha="left", va="top", fontsize=7.5)
    b.text(zV * 2.6, 0.72, r"$87.7\%$ accumulated" + "\n" + "in the screened regime", ha="left", va="center", fontsize=7.5)
    b.text(zs * 1.12, 0.6, "unscreened window", color="#2b7a5a", fontsize=7, ha="center", va="center", rotation=90)
    b.set_xlabel(r"redshift $z$")
    b.set_ylabel(r"cumulative share of $I_1$ above $z_*$")
    panel_label(b, "(b)", x=-0.04, y=1.08)
    # (c) scale dilemma
    kk = np.geomspace(2e-4, 20.0, 200)
    c.loglog(kk, kk * (I1 / (100 * 0.7 / 299792.458)) / np.log(10.0), color=VIOLET, lw=1.9)
    c.axhline(1.0, color=GREY, lw=0.9, ls=(0, (3, 2)))
    c.axvspan(k_g1, k_rec, color="#dbe7f7", alpha=0.9, lw=0)
    c.axvline(k_rec, color=BLUE, lw=0.9); c.axvline(k_g1, color=BLUE, lw=0.9)
    c.set_xlim(2e-4, 20); c.set_ylim(2e-2, 3e4)
    c.text(k_g1 * 0.86, 2e4, r"illustrative cut: $\Gamma\leq10$", ha="right", va="top", fontsize=7, color=INK2, rotation=90)
    c.text(k_rec * 1.18, 2e4, r"horizon entry at $z_{\rm rec}$", ha="left", va="top", fontsize=7, color=INK2, rotation=90)
    c.text(np.sqrt(k_g1 * k_rec), 0.025, "super-horizon", ha="center", va="bottom", fontsize=6.5, color=BLUE, rotation=90)
    kp = [r["k"] for r in km["rows"]]; lp = [r["log10Gamma"] for r in km["rows"]]
    c.plot(kp, lp, "o", color=VIOLET, ms=4.3, mec="white", mew=0.5, zorder=5)
    c.annotate(r"$k=0.05$ Mpc$^{-1}$: $\log_{10}\Gamma=87.87$", xy=(0.05, 87.87), xytext=(0.012, 0.25), fontsize=7, ha="left", va="center", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    c.set_xlabel(r"wavenumber $k$ [Mpc$^{-1}$]")
    c.set_ylabel(r"formal growth exponent $\log_{10}\Gamma(k)$")
    panel_label(c, "(c)", x=-0.04, y=1.08)
    save(fig, "fig6_scalar_sector")


# ---------------------------------------------------------------- Fig. 7
def fig7_selfsimilarity():
    th = 6.1
    t = np.tan(np.radians(th)) ** 2
    s2 = np.sin(np.radians(th)) ** 2
    ss = selfsim.branch_constants(theta_deg=th)
    CK.close("Fig7 u*", ss["ustar"], 30.337301, 1e-6)
    CK.close("Fig7 y*", ss["ystar"], 0.705465, 1e-6)
    CK.close("Fig7 c", ss["c"], 0.582199, 1e-6)
    CK.close("Fig7 kappa", ss["kappa"], 0.374329, 1e-6)
    CK.close("Fig7 du*/deps", selfsim.eps_linearity_slope(), -0.626809, 1e-6)
    hl = selfsim.Hstar_law_check()
    CK.close("Fig7 H(z*) law vs full model at the benchmark (+0.046%)", 100 * hl["rel"], 0.046, 3e-3)
    pl = selfsim.plateau_scan()
    CK.true("Fig7 plateau r = H(z*)/m_FP flat to 9e-7 on m_FP/H0 = 500-8000 (max |r-c| = %.1e)"
            % max(abs(r - ss["c"]) for _, r, _ in pl), max(abs(r - ss["c"]) for _, r, _ in pl) < 1.0e-6)
    asy = selfsim.X_asymptote()
    CK.close("Fig7 X(u=1e8) -> -1", asy["X"], -1.0, 4e-6)

    cm = plt.get_cmap("viridis")
    ms = [6, 12, 25, 50, 100, 200]
    zgrid = np.geomspace(0.01, 3000.0, 360)
    ugrid = np.geomspace(0.3, 1e8, 300)
    Xu = np.array([selfsim.X_of_u(t, 0.0, 3.0, u, 0.0) for u in ugrid])

    fig = plt.figure(figsize=(6.5, 5.3), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, 0]); d = fig.add_subplot(gs[1, 1])
    zst = {}
    for i, mf in enumerate(ms):
        mm = Model(theta_deg=th, mFP_H0=float(mf), Omega_L=0.699, Omega_r=9.08e-5)
        Om = mm.solve_Om0()
        Xz = []
        for z in zgrid:
            try:
                Xz.append(scalarsec.X_of_z(mm, z, Om))
            except Exception:
                Xz.append(np.nan)
        Xz = np.array(Xz)
        col = cm(0.05 + 0.82 * i / (len(ms) - 1))
        u = Om * (1 + zgrid) ** 3 / (s2 * mf ** 2)
        a.plot(zgrid, Xz, color=col, lw=1.4)
        b.plot(u, Xz, color=col, lw=1.4)
        zst[mf] = selfsim.zstar_any(mm, Om)
        a.text(zst[mf], 0.0, "", fontsize=1)
        a.plot([zst[mf]], [0.0], "o", color=col, ms=4, mec="white", mew=0.5, zorder=5)
        a.plot([], [], color=col, lw=1.4, label=r"$m_{\rm FP}/H_0=%d$" % mf)
    for ax_ in (a, b):
        ax_.set_xscale("log"); ax_.axhline(0, color=INK2, lw=0.7)
        ax_.set_ylim(-1.35, 1.15)
    a.set_xlim(0.01, 3000)
    a.set_xlabel(r"redshift $z$"); a.set_ylabel(r"$X(z)$ (full model)")
    a.legend(loc="lower left", ncol=1, fontsize=6.5, title=r"dots: $z_*$", title_fontsize=6.5, borderaxespad=0.3, labelspacing=0.25, handlelength=1.4)
    panel_label(a, "(a)", x=-0.02, y=1.08)
    b.plot(ugrid, Xu, color=INK, lw=1.0, ls=(0, (4, 2)), zorder=1)
    b.axvline(ss["ustar"], color=GREY, lw=0.9, ls=(0, (3, 2)))
    b.axhline(-1.0, color=GREY, lw=0.8, ls=(0, (1, 2)))
    b.set_xlim(0.3, 1e9)
    b.text(ss["ustar"] * 1.3, 0.95, r"$u_*=30.337$", fontsize=7.5, ha="left")

    b.text(3e3, 0.55, "dashed: universal $X(u)$\n($\\epsilon=\\rho=0$)", fontsize=7, color=INK2, ha="left", va="center")
    b.set_xlabel(r"$u=\Omega_{\rm m}(1+z)^3/(\sin^2\theta\,m_{\rm FP}^2)$")
    b.set_ylabel(r"$X$")
    panel_label(b, "(b)", x=-0.02, y=1.08)

    # (c) boundary vs mass: scan + master equation
    S = load("scalar_mfp_scan.json")
    mfs = np.array([r["mFP_H0"] for r in S]); zst_s = np.array([r["zstar"] for r in S])
    zV_s = np.array([r["zV"] for r in S])
    sel = (mfs >= 10)
    slope = np.polyfit(np.log(mfs[sel]), np.log1p(zst_s[sel]), 1)[0]
    CK.close("Fig7 slope of ln(1+z*) vs ln m_FP on m_FP/H0=10-200 (paper 0.6665)", slope, 0.6665, 1.5e-3)
    mgrid = np.geomspace(3.0, 220.0, 22)
    zmast = np.array([selfsim.zstar_master(mFP_H0=float(x)) for x in mgrid])
    c.loglog(mgrid, 1 + zmast, color=ORANGE, lw=1.6)
    c.loglog(mfs, 1 + zst_s, "o", color=BLUE, ms=4.5, mec="white", mew=0.5, zorder=5)
    c.loglog(mfs, 1 + zV_s, "-", color=GREY, lw=0.9, ls=(0, (3, 2)))
    mref = np.array([8.0, 150.0])
    c.loglog(mref, 2.2 * (mref / 8.0) ** (2.0 / 3.0) * 1.0, color=INK2, lw=0.8, ls=(0, (1, 2)))
    from matplotlib.lines import Line2D
    c.legend(handles=[Line2D([], [], marker="o", color=BLUE, ls="", ms=4.5, label="full model"), Line2D([], [], color=ORANGE, lw=1.6, label="master Eq. (13)"), Line2D([], [], color=GREY, lw=0.9, ls=(0, (3, 2)), label=r"$1+z_V$"), Line2D([], [], color=INK2, lw=0.8, ls=(0, (1, 2)), label=r"slope $2/3$")], loc="upper left", fontsize=7)

    c.set_xlim(1.3, 260); c.set_ylim(1.0, 100)
    c.set_xlabel(r"$m_{\rm FP}/H_0$"); c.set_ylabel(r"$1+z_*$")
    c.text(0.97, 0.06, r"measured slope $%.4f$" % slope, transform=c.transAxes, ha="right", fontsize=7.5)
    panel_label(c, "(c)", x=-0.04, y=1.08)
    # (d) H(z*)/m_FP and the law
    Hs = []
    for r in S:
        mm = Model(theta_deg=th, mFP_H0=float(r["mFP_H0"]), Omega_L=0.699, Omega_r=9.08e-5)
        Hs.append(np.sqrt(mm.E2_bim(r["zstar"], mm.solve_Om0())) / r["mFP_H0"])
    Hs = np.array(Hs)
    mg = np.geomspace(1.3, 9000.0, 200)
    law = np.sqrt(ss["c"] ** 2 + ss["kappa"] * 0.699 / mg ** 2)
    d.semilogx(mg, law, color=ORANGE, lw=1.6)
    d.semilogx(mfs, Hs, "o", color=BLUE, ms=4.5, mec="white", mew=0.5, zorder=5)
    d.semilogx([p[0] for p in pl], [p[1] for p in pl], "s", mfc="white", mec=VIOLET, ms=4.6, zorder=5)
    d.axhline(ss["c"], color=INK, lw=0.8, ls=(0, (4, 2)))
    d.text(9000, ss["c"] + 0.015, r"$c=0.582199$", ha="right", va="bottom", fontsize=7.5)
    d.legend(handles=[Line2D([], [], color=ORANGE, lw=1.6, label=r"law (11)"), Line2D([], [], marker="o", color=BLUE, ls="", ms=4.5, label=r"full model ($\Omega_r\neq0$)"), Line2D([], [], marker="s", mfc="white", mec=VIOLET, ls="", ms=4.6, label=r"$\Omega_r=0$, deep regime")], loc="upper right", fontsize=7)


    d.set_ylim(0.54, 1.32); d.set_xlim(1.3, 1.2e4)
    d.set_xlabel(r"$m_{\rm FP}/H_0$"); d.set_ylabel(r"$r=H(z_*)/m_{\rm FP}$")
    panel_label(d, "(d)", x=-0.04, y=1.08)
    save(fig, "fig7_selfsimilarity")


FIGS = {1: fig1_background, 2: fig2_mechanism, 3: fig3_anchor_profile,
        4: fig4_anchor_plane, 5: fig5_cmb_transplant, 6: fig6_scalar_sector,
        7: fig7_selfsimilarity}


def main(argv):
    fs.apply_style()
    which = [int(a) for a in argv] or sorted(FIGS)
    for k in which:
        print("figure", k, "...", flush=True)
        FIGS[k]()
    return CK.write(RESDIR / ("figure_checks.txt" if not argv else "figure_checks_partial.txt"))


if __name__ == "__main__":
    sys.exit(1 if main(sys.argv[1:]) else 0)
