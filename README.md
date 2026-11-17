# Positron cooling and annihilation in noble gases

How does a swarm of positrons cool when its slowest particles are also the most likely to annihilate? Connor Louw’s MSci dissertation at Queen’s University Belfast (2026), supervised by Dr Gleb Gribakin, studies this competition in helium, neon, argon, krypton and xenon using a radial Fokker–Planck equation. The Python scripts reproduce the xenon example, where preferential loss of slow positrons changes the momentum distribution of the survivors.

[Read the dissertation (PDF)](Fokker_Planck_Dissertation.pdf) · [Original project brief (PDF)](PositronCoolingProject.tex.pdf)

## Run the xenon example

Download and extract the repository. The example uses Python 3.12 with NumPy and SciPy. Both input tables are built from values included in the scripts.

**In Spyder**, open the files in `scripts/` and run each with **F5**, in this order.

1. [generate_xe_sigma_m.py](scripts/generate_xe_sigma_m.py) builds the momentum-transfer cross section from the xenon scattering phase shifts.
2. [generate_xe_zeff.py](scripts/generate_xe_zeff.py) builds the momentum-dependent annihilation parameter, $Z_{\mathrm{eff}}(p)$.
3. [xenon_cooling.py](scripts/xenon_cooling.py) reads both tables from `data/` and evolves the positron distribution. Once the tables exist, you can rerun this script on its own.

**From a terminal**, run these commands from the repository’s top-level folder.

```bash
python -m pip install numpy scipy
python scripts/generate_xe_sigma_m.py
python scripts/generate_xe_zeff.py
python scripts/xenon_cooling.py
```

## Results

The solver writes `results/xe_diagnostics.txt`. Its five columns give the time–density product $\tau$ (ns amagat), surviving fraction, mean annihilation rate (per ns amagat), mean $Z_{\mathrm{eff}}$, and mean kinetic energy (Hartree). The three means are taken over the surviving positrons.

For the default 293 K xenon run, the final row gives a **surviving fraction of $1.638\times10^{-13}$** and **mean $Z_{\mathrm{eff}}$ of 397.12**. The mean energy stays above the cooling threshold of $1.01\times(3k_{\mathrm B}T/2)$. The original solver records diagnostics after each sampled update using the time at the start of that step. The row labelled 500 therefore represents 500.001 ns amagat.

The scattering and annihilation inputs come from [Green, Ludlow and Gribakin, *Physical Review A* **90**, 032712 (2014)](https://doi.org/10.1103/PhysRevA.90.032712).
