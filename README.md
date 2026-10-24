# Positron cooling in noble gases

Companion repository for **Connor Louw's MSci project**, *Fokker-Planck
Equation for Positron Cooling and Annihilation*, Queen's University Belfast
(2026), supervised by Dr Gleb Gribakin.

The thesis studies positron cooling and annihilation in He, Ne, Ar, Kr and Xe
using a radial Fokker-Planck equation with many-body-theory scattering inputs.
The xenon study examines how preferential annihilation of slow positrons
affects the distribution of the surviving population.

**[Read the dissertation](Fokker_Planck_Dissertation.pdf).** The copy archived
here is dated **17 April 2026** on its title page. This repository contains the
report and a runnable xenon example; it is a partial code archive, not a
complete reproduction package for every thesis figure or its eigenanalysis.
The later molecular-gas work is a separate project.

## Contents

| File | Purpose |
| --- | --- |
| [Fokker_Planck_Dissertation.pdf](Fokker_Planck_Dissertation.pdf) | MSci report and full scientific discussion |
| [PositronCoolingProject.tex.pdf](PositronCoolingProject.tex.pdf) | Original project brief |
| [scripts/xenon_cooling.py](scripts/xenon_cooling.py) | Xenon cooling, annihilation and diagnostics |
| [scripts/generate_xe_sigma_m.py](scripts/generate_xe_sigma_m.py) | Build the momentum-transfer cross-section table |
| [scripts/generate_xe_zeff.py](scripts/generate_xe_zeff.py) | Build the annihilation parameter table |
| [CITATION.cff](CITATION.cff) | Citation metadata for GitHub's “Cite this repository” button |

The generators create `data/`; the solver creates `results/`. These generated
files are ignored by Git and can be recreated with the commands below.

## Run the xenon example

Use **Python 3.12** for the tested environment. Download or clone this repository,
then run these commands from its top-level folder:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/generate_xe_sigma_m.py
python scripts/generate_xe_zeff.py
python scripts/xenon_cooling.py
```

On Windows, activate the environment with `.venv\Scripts\activate` instead.
In Spyder, use an interpreter with the listed dependencies, run the two
generators once, then run `scripts/xenon_cooling.py`. The scripts locate their
input and output folders relative to the repository, independent of Spyder's
working directory. Edit the grid, gas temperature and time settings near the
top of the solver as in the original project.

The default run uses 400 momentum intervals (401 nodes), a cutoff of
0.6 a.u., a gas temperature of 293 K and a timestep of 0.001 ns amagat.
Its requested duration is 500 ns amagat. Both input tables are generated from
the values included in the scripts; no separate data download is needed.
Re-running the scripts replaces their generated tables or diagnostics.

## Output and conventions

The solver prints the surviving fraction and mean annihilation parameter, and
writes `results/xe_diagnostics.txt`. It does not create plots automatically.
With the default settings, the surviving fraction is approximately
`1.63792360e-13` and the final mean `Zeff` is `397.12032549`; the mean-energy
cooling threshold is not reached. The tidied solver's complete diagnostics
file was checked to be identical to the original script's output with the same
input tables.

The columns are:

| Column | Meaning | Units |
| --- | --- | --- |
| `tau` | Original recorded time-density label | ns amagat |
| `N(tau)` | Surviving population, normalised initially to one | dimensionless |
| `lambda_bar(tau)` | Mean annihilation rate of survivors | (ns amagat)⁻¹ |
| `Zbar_eff(tau)` | Mean annihilation parameter of survivors | dimensionless |
| `Ebar(tau)` | Mean kinetic energy of survivors | Hartree |

Momentum is in atomic units and the momentum-transfer cross section is in
Bohr-radius squared. `F` is the radial momentum density: its discrete population
is `sum(F) * dp`. The initial array is proportional to momentum on the finite
grid, with zero endpoint values and one discrete normalisation. The cutoff is
a numerical choice, not the physical positronium-formation threshold.

The original timing convention is retained: `int(tau_end / dt) + 1` updates
are performed, and diagnostics are labelled `step * dt` after the update.
Consequently, a saved row's physical time is `tau + dt`; the default final
state is at 500.001 ns amagat. The printed “cooling threshold” uses mean energy
within 1% of `3*kB*Tgas/2`; it is not an annihilation-shoulder time.

## Numerical scope and provenance

The example retains the thesis-era endpoint grid, centred elastic flux,
forward-Euler transport update and exponential annihilation step. This tidy-up
changes names, documentation and file locations, not the numerical equations,
physical parameters or time convention. It also restores the existing
momentum-transfer generator from the supplied thesis code archive.

The original `xenon code.py` is now `scripts/xenon_cooling.py`, and
`Xenon MBT.py` is now `scripts/generate_xe_zeff.py`. Earlier versions remain
available in Git history. The PDFs are unchanged.

The input values are from D. G. Green, J. A. Ludlow and G. F. Gribakin,
“Positron scattering and annihilation on noble-gas atoms,” *Physical Review A*
**90**, 032712 (2014), [doi:10.1103/PhysRevA.90.032712](https://doi.org/10.1103/PhysRevA.90.032712).
The cross-section generator uses the tabulated phase shifts, a scattering-length
anchor and the higher-partial-wave polarisation tail. The annihilation generator
adds a fitted zero-momentum anchor and samples a PCHIP interpolant. The solver
linearly interpolates the saved tables.

Reproducing this example does not establish convergence for every observable
or a comparison with fully matched Monte Carlo inputs. The later CF₄ solver's
convergence study does not apply to this older atomic discretisation. Use the
thesis for the broader results and their interpretation.

## Citing the thesis

For the earlier noble-gas study in the PRA or other work, cite the report:

> C. Louw, *Fokker-Planck Equation for Positron Cooling and Annihilation*,
> MSci project report, Queen's University Belfast (2026).

```bibtex
@mastersthesis{Louw2026PositronCooling,
  author = {Louw, Connor},
  title = {Fokker-Planck Equation for Positron Cooling and Annihilation},
  school = {Queen's University Belfast},
  type = {MSci project report},
  year = {2026},
  url = {https://github.com/clouw2003/Fokker-Planck-for-Positron-cooling-and-annihilation/blob/main/Fokker_Planck_Dissertation.pdf}
}
```

If using the code, also record the exact Git commit used. Cite the Green,
Ludlow and Gribakin paper for its scattering and annihilation inputs.
