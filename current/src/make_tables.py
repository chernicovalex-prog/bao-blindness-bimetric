"""Generate all paper tables from the pre1 pipeline into results/.
Every number in the paper comes from this script (or checks.py gates).
"""
import numpy as np
import json, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(SRC / "probes"))

from core import Model
import scalarsec, bao, sne, cmb

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)


def wtsv(name, header, rows):
    with open(OUT / name, "w") as f:
        f.write("# " + " | ".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    print("wrote", name)


def main():
    m = Model()
    Om_bim = m.solve_Om0()
    B = m.B_params()

    # --- scalar sector ---
    s = scalarsec.run_scalar()
    gamma_rows = []
    for (k, h), lg in sorted(s["log10_gamma"].items()):
        gamma_rows.append(["%.3f" % k, "%.3f" % h, "%.4f" % lg])
    wtsv("table_gamma.tsv", ["k [Mpc^-1]", "h", "log10 Gamma"], gamma_rows)

    # --- BAO ---
    bao.load_desi()
    bl = bao.baseline()
    ab = bao.scan_alpha_beta()
    mfp = bao.scan_mFP()
    th = bao.scan_theta()
    bao_rows = [
        ["baseline", "%.5f" % bl["dchi2"], "%.5f" % bl["Om_bim"], "%.5f" % bl["Om_lcdm"]],
        ["(alpha,beta) min", "%.5f" % min(d for d, _ in ab.values() if isinstance(d, float)),
         "-", "-"],
        ["(alpha,beta) max", "%.5f" % max(d for d, _ in ab.values() if isinstance(d, float)),
         "-", "-"],
        ["degenerate (-1,1)", "NO_FINITE_BRANCH", "-", "-"],
    ]
    for mv, (d, om) in mfp.items():
        bao_rows.append(["m_FP=%g" % mv, "%.5f" % d, "%.5f" % om, "-"])
    for t, (dsc, dfx, om) in th.items():
        bao_rows.append(["theta=%.1f self-consistent" % t, "%.5f" % dsc, "%.5f" % om, "-"])
        bao_rows.append(["theta=%.1f fixed-Om" % t, "%.5f" % dfx, "-", "-"])
    wtsv("table_bao.tsv", ["scan", "Delta chi^2", "Om_bim", "Om_lcdm"], bao_rows)

    # --- SNe ---
    sne.load_pantheon()
    rs = sne.run(0.01)
    rs23 = sne.run(0.023)
    sig = sne.sigma_int_scan()
    sne_rows = [
        ["z>0.01 diag", rs["N"], "%.4f" % rs["diag"]["Om_lcdm"],
         "%.3f" % rs["diag"]["chi2_min"], "%.3f" % rs["diag"]["dchi2"],
         "%.3f" % rs["diag"]["pure"], "%.3f" % rs["diag"]["chi2_dof"]],
        ["z>0.01 full", rs["N"], "%.4f" % rs["full"]["Om_lcdm"],
         "%.3f" % rs["full"]["chi2_min"], "%.3f" % rs["full"]["dchi2"],
         "%.3f" % rs["full"]["pure"], "%.3f" % rs["full"]["chi2_dof"]],
        ["z>0.023 full", rs23["N"], "%.4f" % rs23["full"]["Om_lcdm"],
         "%.3f" % rs23["full"]["chi2_min"], "%.3f" % rs23["full"]["dchi2"],
         "%.3f" % rs23["full"]["pure"], "%.3f" % rs23["full"]["chi2_dof"]],
    ]
    for sv, d in sig.items():
        sne_rows.append(["sigma_int=%.2f diag" % sv, rs["N"], "%.4f" % d["Om_lcdm"],
                         "%.3f" % d["chi2_min"], "%.3f" % d["dchi2"],
                         "%.3f" % d["pure"], "%.3f" % d["chi2_dof"]])
    wtsv("table_sne.tsv",
         ["variant", "N", "Om_lcdm", "chi2_min", "Delta chi2", "pure", "chi2/dof"],
         sne_rows)

    # --- CMB ---
    chi2_ref, p_ref = cmb.fit_lcdm_reference()
    _, _, scan = cmb.theta_scan()
    cmb_rows = []
    for t, d in scan.items():
        cmb_rows.append(["%.1f" % t, "%.5f" % d["chi2"], "%.5f" % d["dchi2"],
                         "%.5f" % d["h"], "%.5f" % d["Om"], "%.5f" % d["R"],
                         "%.4f" % d["la"]])
    wtsv("table_cmb_theta_scan.tsv",
         ["theta", "chi2_bim", "Delta chi2", "h", "Om", "R", "l_a"], cmb_rows)
    split = {str(t): cmb.split_test(t) for t in (10.0, 30.0)}
    wtsv("table_cmb_split.tsv", ["theta", "full", "zero_high", "zero_low"],
         [[t, "%.5f" % v["full"], "%.5f" % v["zero_high"], "%.5f" % v["zero_low"]]
          for t, v in split.items()])

    summary = dict(
        Om_bim=Om_bim, B=B.tolist(), zstar=s["zstar"], Hzstar=s["Hzstar"],
        I1=s["I1"], X={str(k): v for k, v in s["X"].items()},
        gamma={"%g/%g" % k: v for k, v in s["log10_gamma"].items()},
        bao=bl, bao_theta={str(k): v for k, v in th.items()},
        sne=rs, sne_023=rs23,
        cmb_ref_chi2=chi2_ref, cmb_ref_params=p_ref.tolist(),
        cmb_scan={str(k): v for k, v in scan.items()}, cmb_split=split,
    )
    with open(OUT / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("wrote summary.json")


if __name__ == "__main__":
    main()
