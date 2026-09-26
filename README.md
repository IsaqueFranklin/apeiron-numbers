# Apeiron: vacancy kinetic Monte Carlo on a chemically disordered FCC lattice

This repository simulates one vacancy moving through an approximately equiatomic Co, Ni, Cr, Fe, Mn alloy. Its main research output is the physical waiting-time series, `dt_1, dt_2, ...`, produced by rejection-free kinetic Monte Carlo (KMC). Every event also records all twelve possible transitions, so the energy and rate calculations can be audited.

The energetic model is motivated by the vacancy migration-barrier statistics reported by [Thomas and Patala (2020)](https://doi.org/10.1016/j.actamat.2020.06.022). Their 2,971 NEB calculations yielded an approximately Gaussian migration-barrier distribution with mean 0.81 eV and standard deviation 0.32 eV. Here, barriers are derived from a well energy at the starting vacancy state and one symmetric saddle energy for the transition:

```text
E_ij = s_ij - w_i
r_ij = nu exp(-E_ij / (k_B T))
```

The independent energy scales are `w ~ N(0, 0.110^2)` eV and `s ~ N(0.81, 0.304^2)` eV. Saddles below either endpoint well are rejected and redrawn. This slightly truncates the model distribution, so the realized moments should be inspected rather than assumed equal to the untruncated values.

## Install

Python 3.10 or newer is recommended. No ASE or AFLOW installation is needed.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Run

```bash
python main.py --jumps 10000 --temperature 1273 --cells 8 \
  --attempt-frequency 1e13 --seed 2020 --output results/run-1273K --plot
```

All options are available through `python main.py --help`. `--cells 8` creates `4 × 8^3 = 2048` FCC sites including one vacancy. At least three cells per axis are required. `--seed` reproduces lattice composition, energy mapping, and KMC draws. Specify a different output directory for each run. `--temperature` is in kelvin, `--attempt-frequency` is in inverse seconds, and `--jumps` is the number of vacancy swaps.

Before each run, the program draws one million independent model barriers by default. It prints their mean and standard deviation and stores the full summary. This mandatory preflight can be shortened for quick development checks with `--preflight-samples 10000`; use the default count for research runs. The preflight checks the ensemble model, while the `realized_offered_barriers` summary describes the barriers encountered along the actual vacancy path.

Output files:

| File | Contents |
| --- | --- |
| `waiting_times.csv` | Event number and `dt_s`, the primary time series |
| `events.csv` | Origin, destination, moving atom, vacancy location, source well, twelve saddles, twelve barriers, twelve rates, total rate, waiting time, cumulative time, and random draws |
| `summary.json` | Configuration, preflight and realized barrier moments, rejection counts, time, and bit-test results |
| `bits.txt` | Exploratory one-bit-per-event digitization of waiting times |
| `diagnostics.png` | Barrier histograms and waiting-time trace, if `--plot` was set |

The twelve array entries in an event follow `apeiron.model.OFFSETS`; `selected_neighbor` is their zero-based index. `vacancy` is the position after the recorded swap. Coordinates are integer FCC coordinates, with periodic wraparound.

## Statistical checks and interpretation

The included bit analysis compares each waiting time against the median of the same run, then applies a small subset of NIST SP 800-22 Rev. 1a: frequency, block frequency, and runs. The median is estimated from the complete series, so the bits are an exploratory diagnostic with a built-in global balance constraint. The code does not implement the complete NIST suite, and passing its tests would not establish cryptographic security or a physical true random number generator. This simulation uses a seeded pseudorandom generator and its output is reproducible. The central result for the TCC is the waiting-time series and its physical/statistical diagnostics.

Run the automated checks with:

```bash
python -m unittest discover -s tests -v
```

For a detailed explanation of every command option, the KMC equations, and bit generation, see [docs/theory-and-run-guide.md](docs/theory-and-run-guide.md). For the code map and modeling decisions, see [docs/architecture.md](docs/architecture.md).
