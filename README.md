# Apeiron

*Vacancy diffusion and waiting-time simulation in a chemically disordered FCC alloy*

Apeiron follows one vacancy as it exchanges places with atoms in an approximately equiatomic Co, Ni, Cr, Fe, Mn lattice. The simulation uses rejection-free kinetic Monte Carlo (KMC). Its primary output is a time series of physical waiting times, one value for each vacancy jump.

The model is informed by [Thomas and Patala (2020)](https://doi.org/10.1016/j.actamat.2020.06.022), who measured an approximately Gaussian distribution of migration barriers in this alloy. Apeiron represents each barrier through a local vacancy well and a symmetric transition saddle. It records every candidate jump so that the energy and rate calculations can be checked afterward.

## Quick start

Python 3.10 or newer is recommended. Install the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Run 10,000 jumps at 1273 K and write a plot:

```bash
python main.py --jumps 10000 --temperature 1273 --cells 8 \
  --attempt-frequency 1e13 --seed 2020 --output results/run-1273K --plot
```

Use `python main.py --help` to list all options. Each run can have its own temperature, jump count, lattice size, attempt frequency, seed, and output directory. The [run guide](docs/theory-and-run-guide.md) explains every argument in the example.

## What happens in a run

An FCC site has 12 nearest neighbors. At each event, Apeiron calculates a migration barrier and rate for all 12 possible jumps:

```text
E_ij = s_ij - w_i
r_ij = nu exp[-E_ij / (k_B T)]
```

Here `w_i` depends on the vacancy's local chemical environment. The transition saddle `s_ij` is shared by the forward and reverse jump. The vacancy chooses a neighbor with probability proportional to its rate, swaps with that atom, and advances physical time by `ln(1/u) / R`, where `R` is the sum of the 12 rates.

Before the walk, a default check samples one million model barriers and reports their mean, spread, and saddle rejection rate. New saddles below either endpoint well are redrawn to avoid negative barriers. The check can be shortened for development with `--preflight-samples 10000`.

## Results

All files go to the directory supplied by `--output`:

| File | Purpose |
| --- | --- |
| `waiting_times.csv` | Main result: event number and waiting time `dt_s` |
| `events.csv` | Full event log, including the 12 saddles, barriers, and rates considered at each jump |
| `summary.json` | Run settings, barrier statistics, rejection rate, total time, and bit-test results |
| `bits.txt` | Exploratory bits obtained from the waiting-time series |
| `diagnostics.png` | Barrier histograms and waiting-time trace when `--plot` is used |

The independent million-barrier check, all barriers offered during the walk, and the barriers actually selected are reported separately. They are different samples and need not have identical distributions.

## Bit analysis and limits

After a run, Apeiron compares each waiting time with the median waiting time of that same run. A time above the median becomes `1`; every other time becomes `0`. It then applies three tests from NIST SP 800-22 Rev. 1a: monobit frequency, block frequency, and runs.

This bit analysis is exploratory. Using the run's own median forces an almost equal number of zeros and ones, so a good monobit result alone says little. The program uses seeded pseudorandom draws; its output is reproducible and is not a physical true random number generator. The waiting-time series is the research result to analyze first.

## Documentation and tests

Read the [theory and run guide](docs/theory-and-run-guide.md) for the derivation of the KMC timing rule, a parameter-by-parameter explanation of the command, and the exact bit conversion. Read the [architecture guide](docs/architecture.md) for the lattice representation, chemical energy keys, module roles, and modeling limitations.

Run the automated tests with:

```bash
python -m unittest discover -s tests -v
```
