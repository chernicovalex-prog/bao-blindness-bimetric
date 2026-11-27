# Re-do the Omega_Lambda profile at theta=6.1 with SNe on the STANDARD Pantheon+ convention (zHD distances).
# Uses the author's BAO/CMB/SNe modules for everything else, so only the SNe redshift convention changes.
import sys, numpy as np, pathlib
sys.path.insert(0,'.'); sys.path.insert(0,'probes')
from core import Model, OMEGA_R_NU
import probes.sne as S
from probes.sne import *
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
S.load_pantheon()
DF=S.DF
mask = DF["zHD"]>0.01
cov=DF["COV_FULL"][np.ix_(mask,mask)]; Ci=np.linalg.inv(cov)
def chi(mu):
    r=DF["m_b_corr"][mask]-mu; Cr=Ci@r; return float(r@Cr-Cr.sum()**2/Ci.sum())
zd=DF["zHD"][mask]; zh=DF["zHEL"][mask]
def mu_of(E2arr):
    ic=cumulative_trapezoid(1/np.sqrt(E2arr),S.Z_FINE,initial=0)
    DL=(1+zh)*S.C_H0*np.interp(zd,S.Z_FINE,ic); return 5*np.log10(DL)+25
from scipy.optimize import minimize_scalar
def lc(Om): return mu_of(Om*(1+S.Z_FINE)**3+OMEGA_R_NU*(1+S.Z_FINE)**4+(1-Om-OMEGA_R_NU))
r=minimize_scalar(lambda Om:chi(lc(Om)),bounds=(0.2,0.5),method='bounded',options={'xatol':1e-10})
print("zHD LCDM: Om=%.4f chi2=%.3f"%(r.x,r.fun))
for OL in (0.64,0.66,0.67,0.68,0.69,0.696,0.699,0.705,0.72):
    m=Model(Omega_L=OL,Omega_r=OMEGA_R_NU); Om=m.solve_Om0()
    E2=np.array([m.E2_bim(z,Om) for z in S.Z_FINE])
    print("OL=%.3f Om=%.4f  SNe(zHD) dchi2 = %.3f"%(OL,Om,chi(mu_of(E2))-r.fun))
