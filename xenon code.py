# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 15:30:47 2026

@author: clouw
"""

import numpy as np

# ============================================================
# Xenon cooling and annihilation solver
# ------------------------------------------------------------
# This script evolves the radial momentum distribution F(p, tau)
# for positrons cooling in xenon, including an annihilation sink.
#
# Input files:
#   - xe_sigma_m_01res.txt
#   - xe_zeff_total_on_grid.txt
#
# Time variable:
#   tau = t(ns) * n_g(amagat)
# ============================================================


# ----------------------------
# 1. Momentum grid
# ----------------------------
pmax = 0.6
Np = 400

p = np.linspace(0.0, pmax, Np + 1)
dp = p[1] - p[0]


# ----------------------------
# 2. Physical parameters
# ----------------------------
Tgas = 293.0
kT_au = 3.1668e-6 * Tgas
D0 = 5.21e-4

# Xenon mass in electron masses
M_Xe = 131.293 * 1822.888486

# Diffusion prefactor appearing in the FP equation
alpha = D0 * (Tgas / M_Xe)


# ----------------------------
# 3. Load transport and annihilation data
# ----------------------------
# Momentum-transfer cross section sigma_m(p)
sigma_data = np.loadtxt("xe_sigma_m_01res.txt")
p_sigma = sigma_data[:, 0]
sigma_tab = sigma_data[:, 1]
sigma_node = np.interp(p, p_sigma, sigma_tab)

# Zeff(p)
zeff_data = np.loadtxt("xe_zeff_total_on_grid.txt")
p_zeff = zeff_data[:, 0]
zeff_tab = zeff_data[:, 1]
zeff_node = np.interp(p, p_zeff, zeff_tab)

# Annihilation rate lambda(p) in (ns amagat)^(-1)
r0 = 2.8179403262e-15
c = 2.99792458e8
n_amg = 2.6867805e25

lambda_prefactor = np.pi * r0**2 * c * n_amg * 1.0e-9
lam = lambda_prefactor * zeff_node


# ----------------------------
# 4. Initial condition
# ----------------------------
# Uniform in energy corresponds to F(p,0) proportional to p
F = 2.0 * p / pmax**2
F[0] = 0.0
F[-1] = 0.0

# Normalise to unit total population
F /= np.sum(F) * dp
F0 = F.copy()


# ----------------------------
# 5. Time-stepping parameters
# ----------------------------
tau_end = 500.0
dt = 1.0e-3
nsteps = int(tau_end / dt) + 1


# ----------------------------
# 6. Diagnostics
# ----------------------------
tau_hist = []
N_hist = []
lambda_avg_hist = []
zeff_avg_hist = []
Ebar_hist = []

# Thermal average energy and threshold used for "cooling time"
E_thermal = 1.5 * kT_au
E_thresh = 1.01 * E_thermal

print(f"lambda prefactor = {lambda_prefactor:.6e} (ns amagat)^(-1)")
print(f"lambda range     = [{lam.min():.6e}, {lam.max():.6e}]")
print(f"max(lambda dt)   = {np.max(lam * dt):.6e}")


# ----------------------------
# 7. Main time-stepping loop
# ----------------------------
for step in range(nsteps):

    # Cell-face values
    p_face = 0.5 * (p[:-1] + p[1:])
    F_face = 0.5 * (F[:-1] + F[1:])
    dFdp_face = (F[1:] - F[:-1]) / dp

    # Interpolate sigma_m onto faces
    sigma_face = np.interp(p_face, p, sigma_node)

    # Coefficients in the flux
    G_face = sigma_face * p_face
    A_face = (p_face / kT_au) - (2.0 / p_face)

    # Conservative flux for the cooling operator
    J_face = -alpha * G_face * (dFdp_face + A_face * F_face)

    # No-flux boundary conditions
    J_face[0] = 0.0
    J_face[-1] = 0.0

    # Conservative update for the transport part
    rhs = np.zeros_like(F)
    rhs[1:-1] = -(J_face[1:] - J_face[:-1]) / dp

    F_star = F + dt * rhs
    F_star[0] = 0.0
    F_star[-1] = 0.0

    # Annihilation step
    F = F_star * np.exp(-lam * dt)
    F[0] = 0.0
    F[-1] = 0.0

    # Save diagnostics every 100 steps
    if step % 100 == 0:
        tau = step * dt
        N = np.sum(F) * dp

        if N > 0.0:
            lambda_avg = np.sum(lam * F) * dp / N
            zeff_avg = np.sum(zeff_node * F) * dp / N
            Ebar = np.sum(0.5 * p**2 * F) * dp / N
        else:
            lambda_avg = np.nan
            zeff_avg = np.nan
            Ebar = np.nan

        tau_hist.append(tau)
        N_hist.append(N)
        lambda_avg_hist.append(lambda_avg)
        zeff_avg_hist.append(zeff_avg)
        Ebar_hist.append(Ebar)


# Convert to arrays
tau_hist = np.array(tau_hist)
N_hist = np.array(N_hist)
lambda_avg_hist = np.array(lambda_avg_hist)
zeff_avg_hist = np.array(zeff_avg_hist)
Ebar_hist = np.array(Ebar_hist)


# ----------------------------
# 8. Cooling time from mean energy
# ----------------------------
cooling_time = None

valid = np.isfinite(Ebar_hist)
tau_valid = tau_hist[valid]
Ebar_valid = Ebar_hist[valid]

crossings = np.where(Ebar_valid <= E_thresh)[0]

if len(crossings) > 0:
    i = crossings[0]

    if i == 0:
        cooling_time = tau_valid[0]
    else:
        t1, t2 = tau_valid[i - 1], tau_valid[i]
        E1, E2 = Ebar_valid[i - 1], Ebar_valid[i]

        if E2 != E1:
            cooling_time = t1 + (E_thresh - E1) * (t2 - t1) / (E2 - E1)
        else:
            cooling_time = t2

    print(f"\nCooling threshold reached at tau = {cooling_time:.6f} ns amagat")
else:
    print("\nCooling threshold was not reached in the simulated time window.")


# ----------------------------
# 9. Late-time annihilation diagnostics
# ----------------------------
n_tail = max(5, len(zeff_avg_hist) // 10)

zbar_final = zeff_avg_hist[-1]

print("\nLate-time diagnostics:")
print(f"Final surviving fraction N(tau_end) = {N_hist[-1]:.8e}")
print(f"Final Zbar_eff                     = {zbar_final:.8f}")


# ----------------------------
# 10. Save output for plotting
# ----------------------------
np.savetxt(
    "xe_diagnostics.txt",
    np.column_stack((tau_hist, N_hist, lambda_avg_hist, zeff_avg_hist, Ebar_hist)),
    header="tau  N(tau)  lambda_bar(tau)  Zbar_eff(tau)  Ebar(tau)"
)
