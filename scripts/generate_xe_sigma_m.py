# -*- coding: utf-8 -*-
"""Generate the xenon momentum-transfer table in data/.

Phase shifts and constants: Green, Ludlow and Gribakin, Phys. Rev. A 90,
032712 (2014), Tables I and VII, https://doi.org/10.1103/PhysRevA.90.032712.
Original filename in the supplied thesis archive: MBT data for Xenon.py."""

from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# --- 1. Xenon constants ---
alpha_d = 27.04   # a.u.
a_xe    = -84.5   # a.u. (Table I)

# k-grid: include k=0 anchor, then MBT-supported grid from 0.02 to 0.60
k_grid = np.concatenate(([0.0], np.arange(0.02, 0.61, 0.01)))

# --- 2. MBT Table VII data (Xe) for l=0,1,2 ---
k_mbt  = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60])
d0_mbt = np.array([0.9849, 1.151, 1.130, 1.085, 1.027, 0.7199, 0.4256, 0.1551, -0.0922, -0.3192])
d1_mbt = np.array([2.303e-3, 9.418e-3, 2.178e-2, 4.001e-2, 6.450e-2, 0.2666, 0.4784, 0.5624, 0.5354, 0.4503])
d2_mbt = np.array([3.223e-4, 1.285e-3, 2.856e-3, 5.015e-3, 7.809e-3, 3.572e-2, 9.425e-2, 0.1904, 0.3118, 0.4298])

cs_d0 = CubicSpline(k_mbt, d0_mbt)
cs_d1 = CubicSpline(k_mbt, d1_mbt)
cs_d2 = CubicSpline(k_mbt, d2_mbt)

kmin = float(k_mbt.min())

# --- 3. Build sigma_m(k) ---
sig_m_res = []

for k in k_grid:
    if k == 0.0:
        sig_m_res.append(4.0 * np.pi * a_xe**2)  # threshold anchor
        continue

    kk = max(k, kmin)  # safety clamp (shouldn't be needed since k>=0.02)
    d0 = float(cs_d0(kk))
    d1 = float(cs_d1(kk))
    d2 = float(cs_d2(kk))

    deltas = [d0, d1, d2]

    # polarization tail for l = 3..10
    for l in range(3, 11):
        dl = (np.pi * alpha_d * k**2) / ((2*l - 1) * (2*l + 1) * (2*l + 3))
        deltas.append(dl)

    # momentum-transfer cross section
    s_m = (4.0 * np.pi / k**2) * sum(
        (l + 1) * np.sin(deltas[l] - deltas[l + 1])**2
        for l in range(len(deltas) - 1)
    )
    sig_m_res.append(s_m)

np.savetxt(DATA_DIR / "xe_sigma_m_01res.txt",
           np.column_stack((k_grid, sig_m_res)),
           header="k sig_m")
print(f"Saved {DATA_DIR / 'xe_sigma_m_01res.txt'}")
