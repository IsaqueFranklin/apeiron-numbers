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

For long runs, the terminal reports progress every 10,000 jumps, including elapsed time and an estimated time remaining. Change this with `--progress-interval N`, or set it to `0` to hide jump progress.

To create plots from a run that has already finished, without repeating the simulation:

```bash
python -m apeiron.plotting results/run-1273K
```

This writes a high-resolution `diagnostics.png` and a vector `diagnostics.pdf` in that run's directory.

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
| `diagnostics.png`, `diagnostics.pdf` | Six-panel figure when `--plot` is used or `python -m apeiron.plotting` is run later |

The independent million-barrier check, all barriers offered during the walk, and the barriers actually selected are reported separately. They are different samples and need not have identical distributions.

The diagnostic figure shows offered versus selected migration barriers, a waiting-time histogram and trace, consecutive equal-bit run lengths, the fraction of ones in the blocks used by the test, and short-lag bit autocorrelation. The bit panels include a fair-coin reference for visual comparison; they are exploratory plots, not additional NIST tests.

## Bit analysis and limits

After a run, Apeiron compares each waiting time with the median waiting time of that same run. A time above the median becomes `1`; every other time becomes `0`. It then applies three tests from NIST SP 800-22 Rev. 1a: monobit frequency, block frequency, and runs.

This bit analysis is exploratory. Using the run's own median forces an almost equal number of zeros and ones, so a good monobit result alone says little. The program uses seeded pseudorandom draws; its output is reproducible and is not a physical true random number generator. The waiting-time series is the research result to analyze first.

## Documentation and tests

Read the [theory and run guide](docs/theory-and-run-guide.md) for the derivation of the KMC timing rule, a parameter-by-parameter explanation of the command, and the exact bit conversion. Read the [architecture guide](docs/architecture.md) for the lattice representation, chemical energy keys, module roles, and modeling limitations.

Run the automated tests with:

```bash
python -m unittest discover -s tests -v
```

## Three TCC scenarios

The CLI now accepts multiple scenarios as `NAME:SIGMA_S:SIGMA_W`, in that order:

```bash
python main.py --scenarios H0:0:0 H1:0.152:0.055 H2:0.304:0.110 \
  --jumps 1000000 --temperature 1273 --cells 8 --attempt-frequency 1e13 \
  --seed 2020 --preflight-samples 1000000 --output results/tcc --plot
```

Each scenario receives a fresh alloy and independent energy caches. The same composition seed gives the same initial chemical state; KMC trajectories evolve naturally with their rates. `--sigma-s` and `--sigma-w` configure a single run. `--mu-s` and `--mu-w` default to 0.81 and 0.0 eV. A zero sigma yields the exact mean. Invalid constant-saddle/unbounded-well combinations are rejected.

Outputs include `manifest.json`, source snapshots and SHA-256 hashes, per-scenario summaries, raw independent preflight samples (`preflight_barriers.npy`), incremental event/time CSVs, and bits. `comparacao/metricas.csv` and `resumo.json` distinguish independent, offered and selected barriers. With `--plot`, four comparative figures are saved in black on white as 350 dpi PNG and SVG; individual diagnostics also retain PDF output. Failed or interrupted runs are marked `incomplete` and excluded from completed-run comparisons. Existing run files are never overwritten.

Storage is checked before starting; allow several GB for three million full records. Offered moments/histograms are accumulated in batches of 4,096 events, without a twelve-million-float Python list. The chemical keys are stored as compact exact byte signatures; their original tuple representations are restored before seeded energy draws, preserving the reference H2 trajectory.

The three tests are frequency, block frequency, and runs at alpha=0.01. The block size defaults to max(20, ceil(n/50)), giving 50 full blocks of 20,000 bits for one million bits and no discarded bits. This meets the size recommendations in SP 800-22 section 2.2.7. Use --block-size to specify an explicit valid size. Applicability, statistics, p-values, ties and residual counts are recorded. See [Portuguese implementation guide](docs/funcionamento-codigo.md).

## Local thesis files

The thesis and completed datasets live outside this code repository in the sibling directory `../tcc-apeiron-words`:

- `gerados/experimentos-tcc/`: final Typst sources, PDFs, tables and figures.
- `results/tcc/`: complete H0, H1 and H2 datasets, summaries and audit records.

Code, tests, scripts and documentation remain in this repository. Generated outputs and Python environments are ignored by Git. Use an external directory with `--output` to keep future datasets outside the repository as well.

To audit the existing datasets from this repository:

```bash
python scripts/verify_tcc_results.py ../tcc-apeiron-words/results/tcc
```

`scripts/build_thesis_bundle.py` reads thesis sources from `../tcc-apeiron-words` by default. Its `--thesis-dir` and `--output` options allow other locations. This local thesis directory is not required to run the simulator or its tests.

`scripts/compile_thesis.py` is an optional document helper and additionally requires the Python `typst` package. The simulation dependencies remain those listed in `requirements.txt`.
