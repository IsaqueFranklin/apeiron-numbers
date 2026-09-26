# Theory and run guide

This guide explains the physical model, the exact command-line settings, and how a waiting-time series becomes the exploratory bit sequence. The implementation details and module map are in [architecture.md](architecture.md).

## Reading a run command

```bash
python main.py --jumps 10000 --temperature 1273 --cells 8 \
  --attempt-frequency 1e13 --seed 2020 --output results/run-1273K --plot
```

The backslash at the end of the first line tells a Unix shell to continue the same command on the next line. It is not a simulation parameter. The equivalent command can be written on one line.

| Part | Meaning | Consequence in this example |
| --- | --- | --- |
| `python main.py` | Start the Python command-line program | Construct a lattice, run the model check, simulate, and write files |
| `--jumps 10000` | Number of executed vacancy-atom exchanges | Exactly 10,000 KMC events, waiting times, and raw exploratory bits |
| `--temperature 1273` | Absolute temperature in kelvin | Sets the Boltzmann factor through `k_B T`; changing it changes rates and often the chosen route |
| `--cells 8` | Number of conventional FCC cells along each axis | `4 × 8³ = 2048` sites, consisting of 2,047 atoms and one vacancy |
| `--attempt-frequency 1e13` | Common attempt frequency `nu`, in inverse seconds | `10¹³ s⁻¹` multiplies every transition rate and sets the time scale |
| `--seed 2020` | Integer seed for reproducibility | Fixes the initial shuffled alloy, environment-to-energy mapping, KMC draws, and preflight draws |
| `--output results/run-1273K` | Directory for this run, relative to the repository when invoked there | Saves CSV, JSON, bits, and plot in that directory; a later run using the same path overwrites those files |
| `--plot` | Boolean switch | Write `diagnostics.png`; omit it to skip plotting |

Additional options are shown by `python main.py --help`. `--preflight-samples` defaults to `1000000` independent barrier samples. `--cells` must be at least 3, and the temperature, attempt frequency, preflight count, and number of jumps must be positive. The code has no lattice constant in meters, so positions are dimensionless FCC lattice coordinates. It reports time in seconds, not a physical displacement or diffusion coefficient.

The seed is split into reproducible streams: `seed` initializes the alloy, `seed + 1` assigns chemical energies, `seed + 2` drives KMC choices and waiting times, and `seed + 10` drives the independent barrier preflight. Keeping the seed fixed while changing only `nu` leaves the jump probabilities the same because `nu` is common to all twelve paths; all waiting times scale inversely with `nu`. Changing temperature generally changes both jump probabilities and waiting times.

## Why wells and saddles are separate

For a vacancy state `i`, the local well `w_i` is its energy reference. A transition to neighbor `j` passes through a saddle with energy `s_ij`. The forward and reverse barriers are

```text
E_ij = s_ij - w_i
E_ji = s_ij - w_j
```

The same saddle is used in both directions, while the two endpoint wells can differ. This is why the model cannot simply attach a barrier to each lattice coordinate. A vacancy swap changes which atom occupies each site and may change the local chemical shell even if the vacancy later returns to the same coordinate. In the implementation, a well is a deterministic cached Gaussian draw keyed by the ordered 12-neighbor chemical shell. A saddle is a deterministic cached draw keyed by the moving species and the bridge environment; that key has the same value under reversal.

[Thomas and Patala (2020)](https://doi.org/10.1016/j.actamat.2020.06.022) reported about 2,971 NEB transitions in equiatomic CoCrFeMnNi and an approximately Gaussian migration-barrier distribution with mean 0.81 eV and standard deviation 0.32 eV. Their paired forward and reverse barriers also give standard deviations `sigma_+ = 0.314 eV` for the symmetric component and `sigma_- = 0.078 eV` for the antisymmetric component. If the well and saddle energies are independent, then

```text
E_plus  = (E_ij + E_ji)/2 = s_ij - (w_i + w_j)/2
E_minus = (E_ij - E_ji)/2 = (w_j - w_i)/2
Var(E_plus)  = sigma_s² + sigma_w²/2
Var(E_minus) = sigma_w²/2
```

Thus `sigma_w = sqrt(2) sigma_- ≈ 0.110 eV`, and `sigma_s = sqrt(sigma_+² - sigma_-²) ≈ 0.304 eV`. The code uses `w ~ N(0, 0.110²)` eV and `s ~ N(0.81, 0.304²)` eV. Before caching a new saddle it redraws any candidate below either endpoint well. Both directional barriers are therefore nonnegative. This rejection changes the exact Gaussian moments slightly. It is rejection sampling of the *energy landscape*, not a rejected KMC jump.

The untruncated variance of a forward barrier is `sigma_s² + sigma_w²`, yielding a standard deviation of about 0.323 eV. The required preflight draws one million independent model transitions by default, applies the same nonnegative-barrier rule, and reports its actual moments and rejection fraction. The preflight measures the model ensemble. The barriers offered along one KMC trajectory are a correlated, visited sample; the chosen barriers are additionally weighted toward faster jumps. Those three distributions answer different questions and are reported separately.

## From energy to a KMC event

At a vacancy position `i`, the FCC lattice has twelve nearest neighbors. For each `j`, the code calculates

```text
r_ij = nu exp[-(s_ij - w_i)/(k_B T)]
R    = sum over the twelve r_ij
```

Here `k_B = 8.617333262 × 10⁻⁵ eV/K`; the quotient in the exponential is dimensionless. A higher barrier or lower temperature generally reduces a rate. At a fixed current state, the source well is common to all twelve rates:

```text
r_ij = nu exp[w_i/(k_B T)] exp[-s_ij/(k_B T)]
P(choose j | current state) = r_ij/R
```

The common well factor cancels from these twelve choice probabilities, but it remains in `R` and therefore affects the waiting time. The saddle differences govern which neighbor is favored; the well influences how long the vacancy waits. For a forward/reverse pair with the same attempt frequency, `r_ij/r_ji = exp[(w_i - w_j)/(k_B T)]`, which is the local detailed-balance relation for the two endpoint states.

KMC uses one uniform random value to select a neighbor according to `r_ij/R`. Once selected, the neighboring atom and the vacancy exchange positions. No candidate KMC move is discarded. Another independent uniform value `u` in `(0, 1]` gives the residence time:

```text
P(dt > t | current state) = exp(-R t)
F(t | current state) = 1 - exp(-R t)
dt = -ln(u)/R = ln(1/u)/R
```

The last line is inverse-transform sampling of an exponential waiting time. Its conditional mean is `1/R`. The simulation adds `dt` to cumulative time and repeats from the changed chemical state. This yields `dt_1, dt_2, ..., dt_N`. Because the state and its `R` change and trajectories can revisit local environments, successive waiting times need not be identically distributed or independent when considered without conditioning on the state.

## Exactly how the current bit file is made

After all `N` jumps, the program reads the *entire* waiting-time sequence and computes its sample median `m`. It then writes one bit per event:

```text
b_k = 1 if dt_k > m, otherwise 0
```

For example, if four waiting times in arbitrary units are `[0.2, 0.9, 0.4, 0.1]`, their median is `0.3` and the bits are `0, 1, 1, 0`. `bits.txt` stores these bits without separators. `summary.json` records `m` as `median_dt_s`. The primary scientific observation remains the continuous `waiting_times.csv`; bit digitization is a separate analysis step.

The threshold is estimated from the very same run that is tested. With an even number of distinct waiting times, this construction forces exactly half the bits to be one. Therefore a frequency or monobit test will look perfect by construction, regardless of whether the underlying process is independent. Changing one time can move the median and affect other bits. The median rule also ignores the magnitude of a waiting time once it has been classified. It is useful for an exploratory comparison of temporal patterns, but it does not extract a proven amount of entropy.

The code reports three p-values from a subset of [NIST SP 800-22 Rev. 1a](https://doi.org/10.6028/NIST.SP.800-22r1a): monobit frequency, frequency within 128-bit blocks, and runs. At the configured significance level `alpha = 0.01`, a p-value below `0.01` flags a failed check. The runs test is omitted if its required overall frequency condition fails. The block test is omitted if there are fewer than 128 bits, and streams under 100 bits are marked insufficient. This is not the full NIST battery; one stream and three tests cannot characterize an entropy source. A failed block or runs test may reveal temporal structure even while the median forces the global bit count to be balanced.

Finally, the program calls seeded pseudorandom generators to draw the alloy, energies, jump choices, and waiting times. The resulting file is reproducible from the seed and is **not a physical true random number generator**. Statistical tests can detect some patterns; passing them would not prove unpredictability or cryptographic security. Studying randomness as a separate research question would require a clearly specified extractor, independent validation data, multiple runs, and an entropy model appropriate to the intended application.
