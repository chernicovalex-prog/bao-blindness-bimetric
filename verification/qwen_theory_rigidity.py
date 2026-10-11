import numpy as np

def compute_omega_m(alpha, beta, theta_deg=6.1, m_FP_ratio=25, Omega_Lambda=0.699, Omega_r=5.38e-5):
    theta_rad = np.deg2rad(theta_deg)
    s2 = np.sin(theta_rad)**2
    c2 = np.cos(theta_rad)**2
    m2 = m_FP_ratio**2
    
    # Coefficients of the quartic P(y) at z=0
    # Omega_tot is unknown initially, but we can solve for y using the closure condition iteratively
    # Or simply solve the system:
    # 1. P(y) = 0
    # 2. Omega_m = 1 - Omega_r - Omega_DE(y)
    # 3. Omega_tot = Omega_m + Omega_r
    
    # Substitute Omega_tot = Omega_m + Omega_r into a1
    # a1 = (Omega_m + Omega_r) + Omega_Lambda + ...
    # This makes a1 dependent on Omega_m, which depends on y.
    # It's a coupled system.
    
    # Better approach: Express Omega_m in terms of y from closure, 
    # then substitute into P(y)=0 which contains Omega_tot = Omega_m + Omega_r.
    
    # Omega_DE(y) = Omega_Lambda - s2 * m2 * (1-y) * (1 + alpha*(1-y) + (beta/3)*(1-y)**2)
    # Omega_m(y) = 1 - Omega_r - Omega_DE(y)
    # Omega_tot(y) = Omega_m(y) + Omega_r = 1 - Omega_DE(y)
    
    # Coeffs:
    # a0 = -1/3 * c2 * m2 * (1 + 2*alpha + beta)
    # a1 = Omega_tot(y) + Omega_Lambda + m2 * (c2*(alpha+beta) - s2*(1+alpha+beta/3))
    # a2 = m2 * (-c2*beta + s2*(1+2*alpha+beta))
    # a3 = -(Omega_Lambda + 1/3 * m2 * (c2*(-1+alpha-beta) + 3*s2*(alpha+beta)))
    # a4 = 1/3 * s2 * m2 * beta
    
    # Note: a1 depends on y via Omega_tot(y).
    # P(y) = a0 + a1(y)*y + a2*y**2 + a3*y**3 + a4*y**4 = 0
    
    def residual(y):
        om_de = Omega_Lambda - s2 * m2 * (1-y) * (1 + alpha*(1-y) + (beta/3)*(1-y)**2)
        om_tot = 1 - om_de
        a0 = -1/3 * c2 * m2 * (1 + 2*alpha + beta)
        a1 = om_tot + Omega_Lambda + m2 * (c2*(alpha+beta) - s2*(1+alpha+beta/3))
        a2 = m2 * (-c2*beta + s2*(1+2*alpha+beta))
        a3 = -(Omega_Lambda + 1/3 * m2 * (c2*(-1+alpha-beta) + 3*s2*(alpha+beta)))
        a4 = 1/3 * s2 * m2 * beta
        return a0 + a1*y + a2*y**2 + a3*y**3 + a4*y**4
    
    # Solve for y near 1
    from scipy.optimize import fsolve
    y_sol = fsolve(residual, 1.0)[0]
    
    # Calculate Omega_m
    om_de = Omega_Lambda - s2 * m2 * (1-y_sol) * (1 + alpha*(1-y_sol) + (beta/3)*(1-y_sol)**2)
    om_m = 1 - Omega_r - om_de
    return om_m, y_sol

# Parameters
alphas = [-2, -1, 0, 1, 2]
betas = [1, 2, 3, 4, 5]

print(f"{'Alpha':>5} | {'Beta':>5} | {'Omega_m':>10} | {'Y':>10}")
print("-" * 40)

results = {}
for a in alphas:
    for b in betas:
        om, y = compute_omega_m(a, b)
        results[(a,b)] = om
        if b == 3: # Print row for Beta=3
             print(f"{a:>5} | {b:>5} | {om:.8f} | {y:.8f}")

# Estimate derivatives at alpha=0, beta=3
om_00 = results[(0,3)]
om_p1 = results[(1,3)]
om_m1 = results[(-1,3)]
d_om_da = (om_p1 - om_m1) / 2.0

om_b_p1 = results[(0,4)]
om_b_m1 = results[(0,2)]
d_om_db = (om_b_p1 - om_b_m1) / 2.0

print(f"\nDerivatives at alpha=0, beta=3:")
print(f"d(Omega_m)/d(alpha) ~ {d_om_da:.2e}")
print(f"d(Omega_m)/d(beta)  ~ {d_om_db:.2e}")
