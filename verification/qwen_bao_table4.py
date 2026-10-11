import numpy as np
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d

# ==============================================================================
# ПАРАМЕТРЫ МОДЕЛИ (Högås & Mörtsell 2025, arXiv:2507.03743)
# ==============================================================================
theta_deg = 6.1
theta = np.radians(theta_deg)
m = 25.0
m2 = m**2
alpha = 0.0
beta = 3.0
Omega_Lambda = 0.699
Omega_r = 5.38e-5

cos2 = np.cos(theta)**2
sin2 = np.sin(theta)**2

def get_Omega_DE(y):
    """Плотность темной энергии в биметрической модели."""
    return Omega_Lambda - sin2 * m2 * (1 - y) * (1 + alpha * (1 - y) + (beta / 3.0) * (1 - y)**2)

def get_y(z, Omega_m):
    """Находит физический корень y(z) ∈ (0, 1] квартики A13."""
    Omega_tot_z = Omega_m * (1 + z)**3 + Omega_r * (1 + z)**4
    
    a0 = -(1/3) * cos2 * m2 * (1 + 2*alpha + beta)
    a1 = Omega_tot_z + Omega_Lambda + m2 * (cos2 * (alpha + beta) - sin2 * (1 + alpha + beta/3.0))
    a2 = m2 * (-cos2 * beta + sin2 * (1 + 2*alpha + beta))
    a3 = -(Omega_Lambda + (1/3) * m2 * (cos2 * (-1 + alpha - beta) + 3 * sin2 * (alpha + beta)))
    a4 = (1/3) * sin2 * m2 * beta
    
    def quartic(y):
        return a0 + a1 * y + a2 * y**2 + a3 * y**3 + a4 * y**4
    
    # Физический корень лежит в интервале (0, 1]
    return brentq(quartic, 1e-10, 1.0)

def get_E(z, Omega_m, model='bimetric'):
    """Функция Хаббла E(z) = H(z)/H0."""
    if model == 'lcdm':
        return np.sqrt(Omega_m * (1 + z)**3 + Omega_r * (1 + z)**4 + Omega_Lambda)
    else:
        y = get_y(z, Omega_m)
        Omega_DE = get_Omega_DE(y)
        Omega_tot_z = Omega_m * (1 + z)**3 + Omega_r * (1 + z)**4
        return np.sqrt(Omega_tot_z + Omega_DE)

# ==============================================================================
# ДАННЫЕ DESI DR2 (Table IV, arXiv:2503.14738)
# ==============================================================================
z_data = np.array([0.295, 0.510, 0.510, 0.706, 0.706, 0.934, 0.934, 1.321, 1.321, 1.484, 1.484, 2.330, 2.330])
types = np.array(['DV', 'DM', 'DH', 'DM', 'DH', 'DM', 'DH', 'DM', 'DH', 'DM', 'DH', 'DM', 'DH'])
D_data = np.array([7.942, 13.588, 21.863, 17.351, 19.455, 21.576, 17.641, 27.601, 14.176, 30.512, 12.817, 38.988, 8.632])
err = np.array([0.075, 0.167, 0.425, 0.177, 0.330, 0.152, 0.193, 0.318, 0.221, 0.760, 0.516, 0.531, 0.101])

# Коэффициенты корреляции между DM и DH для одного и того же тресера
corr_vals = np.array([0.0, -0.459, -0.459, -0.404, -0.404, -0.416, -0.416, -0.434, -0.434, -0.500, -0.500, -0.431, -0.431])

# Построение блок-диагональной ковариационной матрицы 13x13
C = np.diag(err**2)
for i in range(1, 13, 2):
    C[i, i+1] = corr_vals[i] * err[i] * err[i+1]
    C[i+1, i] = C[i, i+1]

C_inv = np.linalg.inv(C)

def compute_chi2(Omega_m, model='bimetric'):
    """Вычисляет минимальный χ² и лучший параметр f = c/(H0*rd) через аналитическую маргинализацию."""
    z_grid = np.linspace(0, max(z_data) + 0.1, 1000)
    E_grid = np.array([get_E(z, Omega_m, model) for z in z_grid])
    
    I_DM_grid = np.zeros_like(z_grid)
    I_DM_grid[1:] = cumulative_trapezoid(1.0 / E_grid, z_grid)
    
    I_DM_func = interp1d(z_grid, I_DM_grid, kind='linear', fill_value='extrapolate')
    E_func = interp1d(z_grid, E_grid, kind='linear', fill_value='extrapolate')
    
    I_DM = I_DM_func(z_data)
    E_z = E_func(z_data)
    I_DH = 1.0 / E_z
    
    I = np.zeros_like(z_data)
    for i in range(len(z_data)):
        if types[i] == 'DV':
            I[i] = (z_data[i] * I_DM[i]**2 * I_DH[i])**(1/3)
        elif types[i] == 'DM':
            I[i] = I_DM[i]
        elif types[i] == 'DH':
            I[i] = I_DH[i]
            
    # Аналитическая маргинализация по f: f_best = (I^T C^-1 D) / (I^T C^-1 I)
    I_C_inv = I @ C_inv
    f_best = (I_C_inv @ D_data) / (I_C_inv @ I)
    chi2_min = (D_data @ C_inv @ D_data) - (I_C_inv @ D_data)**2 / (I_C_inv @ I)
    
    return chi2_min, f_best

# ==============================================================================
# РАСЧЁТЫ
# ==============================================================================
# 1. Нахождение Ω_m и y(0) из условия замыкания E²(0) = 1 для биметрической модели
def quartic_residual_y0(y):
    Omega_DE_y = get_Omega_DE(y)
    Omega_m = 1.0 - Omega_r - Omega_DE_y
    Omega_tot_0 = Omega_m + Omega_r
    
    a0 = -(1/3) * cos2 * m2 * (1 + 2*alpha + beta)
    a1 = Omega_tot_0 + Omega_Lambda + m2 * (cos2 * (alpha + beta) - sin2 * (1 + alpha + beta/3.0))
    a2 = m2 * (-cos2 * beta + sin2 * (1 + 2*alpha + beta))
    a3 = -(Omega_Lambda + (1/3) * m2 * (cos2 * (-1 + alpha - beta) + 3 * sin2 * (alpha + beta)))
    a4 = (1/3) * sin2 * m2 * beta
    
    return a0 + a1 * y + a2 * y**2 + a3 * y**3 + a4 * y**4

y0 = brentq(quartic_residual_y0, 1e-10, 1.0)
Omega_m_bim = 1.0 - Omega_r - get_Omega_DE(y0)

print(f"=== БИМЕТРИЧЕСКАЯ МОДЕЛЬ ===")
print(f"y(0) = {y0:.8f}")
print(f"Ω_m (из замыкания) = {Omega_m_bim:.8f}")
chi2_bim, f_bim = compute_chi2(Omega_m_bim, model='bimetric')
print(f"χ²_min = {chi2_bim:.6f}")
print(f"Лучший f = {f_bim:.6f}")

print(f"\n=== ΛCDM КОНТРОЛЬ (Ω_m = 0.2972) ===")
Omega_m_lcdm = 0.2972
chi2_lcdm, f_lcdm = compute_chi2(Omega_m_lcdm, model='lcdm')
print(f"χ²_min = {chi2_lcdm:.6f}")
print(f"Лучший f = {f_lcdm:.6f}")

print(f"\n=== ИТОГОВЫЕ ЧИСЛА ===")
print(f"Δχ² (bimetric − ΛCDM) = {chi2_bim - chi2_lcdm:.6f}")
