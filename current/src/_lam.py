
import sys, json
sys.path.insert(0, "."); sys.path.insert(0, "probes")
import cmb
out = {}
for th in (6.1, 30.0):
    row = []
    for lam in (0.0, 0.25, 0.5, 0.75, 1.0):
        c2, h, wb, Om, vec = cmb.fit_bimetric(th, r_map=lam)
        row.append([lam, float(c2), float(h)])
        print("th=%.1f lam=%.2f chi2=%.4f h=%.5f" % (th, lam, c2, h), flush=True)
    out[str(th)] = row
json.dump(out, open("../results/cmb_lambda_sweep.json", "w"), indent=2)
print("LAMBDA SWEEP DONE")
