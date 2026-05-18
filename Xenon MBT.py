# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 16:04:58 2026

@author: clouw
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

# ------------------------------------------------------------
# Build the xenon total Zeff(p) table used by the solver.
# The p = 0 anchor is obtained from the first three points
# using the virtual-level form Z(p) = K / (x^2 + p^2) + A.
# ------------------------------------------------------------

def three_point_virtual_level_fit(p_vals, z_vals):
    """Exact 3-point solve for Z(p) = K/(x^2 + p^2) + A."""
    p_vals = np.asarray(p_vals, dtype=float)
    z_vals = np.asarray(z_vals, dtype=float)

    p2 = p_vals**2
    p1, p2_, p3 = p2
    z1, z2, z3 = z_vals

    d12 = p1 - p2_
    d23 = p2_ - p3

    num = (z3 - z2) * d12 * z1 - (z2 - z1) * d23 * z3
    den = (z3 - z2) * d12 - (z2 - z1) * d23

    if np.isclose(den, 0.0):
        raise RuntimeError("3-point fit became degenerate.")

    A = num / den

    # Fit 1 / (Z - A) = b + m p^2
    f = 1.0 / (z_vals - A)
    m = (f[0] - f[1]) / (p2[0] - p2[1])
    b = f[0] - m * p2[0]

    K = 1.0 / m
    x2 = b * K
    Z0 = K / x2 + A

    return K, x2, A, Z0


# Tabulated positron momenta
p_tab = np.array([0.02, 0.04, 0.06, 0.08, 0.10,
                  0.20, 0.30, 0.40, 0.50, 0.60], dtype=float)

# s-wave contributions: valence + core
Zs_val = np.array([1227.0, 351.7, 161.0, 93.91, 62.39,
                   18.88, 10.01, 6.618, 4.906, 3.896])

Zs_core = np.array([29.71, 8.289, 3.767, 2.195, 1.463,
                    4.608e-1, 2.623e-1, 1.912e-1, 1.597e-1, 1.453e-1])

# p-wave contributions: valence + core
Zp_val = np.array([1.040e-1, 4.235e-1, 9.763e-1, 1.782, 2.854,
                   11.22, 16.75, 15.84, 13.28, 11.08])

Zp_core = np.array([1.085e-3, 4.428e-3, 1.025e-2, 1.880e-2, 3.032e-2,
                    1.256e-1, 2.023e-1, 2.110e-1, 1.989e-1, 1.896e-1])

# d-wave contribution
Zd = np.array([2.716e-5, 4.350e-4, 2.205e-3, 6.975e-3, 1.704e-2,
               2.713e-1, 1.325, 3.769, 7.384, 10.99])

# Total Zeff(p)
Z_total = (Zs_val + Zs_core) + (Zp_val + Zp_core) + Zd

# Estimate Zeff(0) from the first three points
K, x2, A, Z0 = three_point_virtual_level_fit(p_tab[:3], Z_total[:3])

# Add the p = 0 anchor before interpolation
p_use = np.insert(p_tab, 0, 0.0)
Z_use = np.insert(Z_total, 0, Z0)

# Shape-preserving interpolation on the solver grid
interp = PchipInterpolator(p_use, Z_use)

p_grid = np.arange(0.0, 0.6000001, 0.01)
Z_grid = np.maximum(interp(p_grid), 0.0)

np.savetxt(
    "xe_zeff_total_on_grid.txt",
    np.column_stack((p_grid, Z_grid)),
    header="p  Zeff_total"
)