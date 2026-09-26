"""Command-line simulation, analysis and plotting."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .analysis import barrier_preflight, bits_from_waiting_times, nist_subset
from .model import EnergyLandscape, FCCLattice, KMC


def run(args):
    if args.jumps < 1:
        raise ValueError("jumps must be positive")
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    preflight = barrier_preflight(args.preflight_samples, args.seed + 10)
    print(f"Preflight: {preflight['count']:,} barriers, mean {preflight['mean_eV']:.4f} eV, std {preflight['std_eV']:.4f} eV")
    lattice = FCCLattice(args.cells, args.seed)
    landscape = EnergyLandscape(lattice, args.seed + 1)
    kmc = KMC(lattice, landscape, args.temperature, args.attempt_frequency, args.seed + 2)
    fieldnames = [
        "event", "origin", "destination", "selected_neighbor", "moving_species",
        "vacancy", "source_well_eV", "saddles_eV", "barriers_eV",
        "rates_per_s", "total_rate_per_s", "dt_s", "time_s",
        "selection_u", "waiting_u",
    ]
    dt = []
    offered_barriers = []
    selected_barriers = []
    with (output / "events.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for number in range(1, args.jumps + 1):
            event = kmc.step(number)
            offered_barriers.extend(event["barriers_eV"])
            selected_barriers.append(event["barriers_eV"][event["selected_neighbor"]])
            dt.append(event["dt_s"])
            writer.writerow({
                key: json.dumps(value) if isinstance(value, (list, tuple)) else value
                for key, value in event.items()
            })
    with (output / "waiting_times.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("event", "dt_s"))
        writer.writerows(enumerate(dt, 1))
    bits, threshold = bits_from_waiting_times(dt)
    (output / "bits.txt").write_text("".join(map(str, bits.tolist())) + "\n")
    offered = np.asarray(offered_barriers)
    selected = np.asarray(selected_barriers)
    diagnostics = {
        "parameters": vars(args), "preflight": preflight,
        "realized_offered_barriers": {
            "count": len(offered), "mean_eV": float(offered.mean()),
            "std_eV": float(offered.std()), "min_eV": float(offered.min()),
        },
        "selected_barriers": {
            "count": len(selected), "mean_eV": float(selected.mean()),
            "std_eV": float(selected.std()),
        },
        "landscape": {
            "unique_wells": len(landscape.wells),
            "unique_saddles": len(landscape.saddles),
            "saddle_draws": landscape.saddle_draws,
            "saddle_rejections": landscape.saddle_rejections,
            "rejection_rate": landscape.saddle_rejections / landscape.saddle_draws,
        },
        "total_time_s": kmc.time,
        "digitization": {"rule": "dt > whole-run median", "median_dt_s": threshold},
        "nist_sp800_22_subset": nist_subset(bits),
    }
    (output / "summary.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    if args.plot:
        make_plot(output, offered, selected, np.asarray(dt))
    print(f"Wrote {output / 'waiting_times.csv'}, events.csv, summary.json and bits.txt")


def make_plot(output, offered, selected, dt):
    import os
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".matplotlib-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].hist(offered, bins=60, density=True, alpha=0.65, label="offered barriers")
    axes[0].hist(selected, bins=60, density=True, alpha=0.45, label="selected barriers")
    axes[0].axvline(0.81, color="black", linestyle=":", label="reference mean")
    axes[0].set(xlabel="Migration barrier (eV)", ylabel="Density")
    axes[0].legend()
    axes[1].plot(np.arange(1, len(dt) + 1), dt, linewidth=0.5)
    axes[1].set(xlabel="KMC event", ylabel="Waiting time (s)")
    axes[1].set_yscale("log")
    fig.tight_layout()
    fig.savefig(output / "diagnostics.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Vacancy diffusion on a periodic FCC alloy")
    parser.add_argument("--cells", type=int, default=8, help="conventional FCC cells per axis")
    parser.add_argument("--jumps", type=int, default=10_000)
    parser.add_argument("--temperature", type=float, default=1273.0, help="kelvin")
    parser.add_argument("--attempt-frequency", type=float, default=1e13, help="per second")
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--preflight-samples", type=int, default=1_000_000)
    parser.add_argument("--output", default="results/run")
    parser.add_argument("--plot", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
