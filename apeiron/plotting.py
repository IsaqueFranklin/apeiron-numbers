"""Publication-ready exploratory plots for a completed KMC run."""

import argparse
import csv
import json
import os
from pathlib import Path

import numpy as np


def _barrier_histograms(events_path, edges):
    offered = np.zeros(len(edges) - 1, dtype=np.int64)
    selected = np.zeros_like(offered)
    offered_batch = []
    selected_batch = []
    with events_path.open(newline="") as handle:
        for event in csv.DictReader(handle):
            barriers = json.loads(event["barriers_eV"])
            offered_batch.extend(barriers)
            selected_batch.append(barriers[int(event["selected_neighbor"])])
            if len(selected_batch) >= 10_000:
                offered += np.histogram(offered_batch, bins=edges)[0]
                selected += np.histogram(selected_batch, bins=edges)[0]
                offered_batch.clear()
                selected_batch.clear()
    offered += np.histogram(offered_batch, bins=edges)[0]
    selected += np.histogram(selected_batch, bins=edges)[0]
    return offered, selected


def plot_results(output, offered=None, selected=None, dt=None, bits=None):
    """Write a six-panel PNG and vector PDF from run arrays or saved files."""
    output = Path(output)
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".matplotlib-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if dt is None:
        with (output / "waiting_times.csv").open(newline="") as handle:
            dt = np.fromiter((float(row["dt_s"]) for row in csv.DictReader(handle)), dtype=float)
    else:
        dt = np.asarray(dt, dtype=float)
    if bits is None:
        raw = (output / "bits.txt").read_bytes().strip()
        bits = np.frombuffer(raw, dtype=np.uint8) - 48
    else:
        bits = np.asarray(bits, dtype=np.uint8)
    if len(dt) != len(bits) or not len(dt):
        raise ValueError("waiting times and bits must have the same nonzero length")

    archive = output / "barrier_histograms.npz"
    edges = np.linspace(0, 3.0, 91)
    if archive.exists() and (offered is None or selected is None):
        with np.load(archive) as hist:
            edges = hist["edges_eV"]
            offered_counts = hist["offered_counts"]
            selected_counts = hist["selected_counts"]
    elif offered is None or selected is None:
        offered_counts, selected_counts = _barrier_histograms(output / "events.csv", edges)
    else:
        offered_counts = np.histogram(offered, bins=edges)[0]
        selected_counts = np.histogram(selected, bins=edges)[0]

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(3, 2, figsize=(8.27, 11.69), constrained_layout=True)
    fig.suptitle(f"Apeiron · {output.name} · {len(dt):,} jumps", fontsize=14)

    ax = axes[0, 0]
    width = np.diff(edges)
    ax.bar(edges[:-1], offered_counts / max(offered_counts.sum(), 1) / width,
           width=width, align="edge", color="white", edgecolor="black", alpha=0.55, label="12 offered per jump")
    ax.step(edges[:-1], selected_counts / max(selected_counts.sum(), 1) / width,
            where="post", color="black", linestyle="--", linewidth=1.5, label="selected jump")
    ax.axvline(0.81, color="black", linestyle=":", label="0.81 eV reference")
    ax.set(xlabel="Migration barrier (eV)", ylabel="Density", title="Barrier selection")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    positive = dt[dt > 0]
    if len(positive):
        time_edges = np.geomspace(positive.min(), positive.max() * (1 + 1e-9), 61)
        ax.hist(positive, bins=time_edges, color="white", edgecolor="black")
        ax.set_xscale("log")
    ax.set(xlabel="Waiting time (s)", ylabel="Events", title="Waiting-time distribution")

    ax = axes[1, 0]
    stride = max(1, len(dt) // 20_000)
    ax.plot(np.arange(1, len(dt) + 1, stride), np.where(dt[::stride] > 0, dt[::stride], np.nan), linewidth=0.45, color="black")
    ax.set_yscale("log")
    ax.set(xlabel="Jump number", ylabel="Waiting time (s)", title="Waiting-time trace")
    if stride > 1:
        ax.text(0.98, 0.98, f"every {stride:,}th point", transform=ax.transAxes,
                ha="right", va="top", fontsize=8)

    ax = axes[1, 1]
    change = np.flatnonzero(bits[1:] != bits[:-1]) + 1
    lengths = np.diff(np.r_[0, change, len(bits)])
    k = np.arange(1, 17)
    observed = np.array([(lengths == i).sum() for i in k[:-1]] + [(lengths >= 16).sum()])
    expected = np.r_[len(bits) / 2.0 ** (k[:-1] + 1), len(bits) / 2.0 ** 16]
    ax.bar(k - 0.2, observed, width=0.4, color="white", edgecolor="black", label="observed")
    ax.bar(k + 0.2, expected, width=0.4, color="white", edgecolor="black", hatch="///", label="fair-coin reference")
    ax.set_yscale("log")
    ax.set_xticks(k[::2], [str(i) for i in k[::2]])
    ax.set(xlabel="Consecutive equal bits (16 = 16+)", ylabel="Number of runs",
           title="Bit run lengths")
    ax.legend(fontsize=8)

    ax = axes[2, 0]
    from .analysis import select_block_size
    summary_path = output / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    block_size = (summary.get("nist_sp800_22_subset", {}).get("block_size")
                  or summary.get("parameters", {}).get("block_size")
                  or select_block_size(len(bits)))
    blocks = len(bits) // block_size
    if blocks:
        fractions = bits[:blocks * block_size].reshape(blocks, block_size).mean(axis=1)
        ax.plot(np.arange(1, blocks+1), fractions, color="black", marker="o", markersize=2,
                linewidth=.6, label="observed")
        ax.axhline(.5, color="black", linestyle="--", label="equal proportions")
        ax.legend(fontsize=8)
    ax.set(xlabel="Block index", ylabel="Fraction of 1s",
           title=f"Local bit balance (M={block_size})")

    ax = axes[2, 1]
    centered = bits.astype(float) - bits.mean()
    variance = np.dot(centered, centered)
    max_lag = min(50, len(bits) - 1)
    if variance and max_lag:
        lags = np.arange(1, max_lag + 1)
        correlation = [np.dot(centered[:-lag], centered[lag:]) / variance for lag in lags]
        ax.bar(lags, correlation, width=0.8, color="black")
        bound = 1.96 / np.sqrt(len(bits))
        ax.axhline(bound, color="black", linestyle="--", linewidth=1)
        ax.axhline(-bound, color="black", linestyle="--", linewidth=1,
                   label="approx. fair-coin 95% bounds")
        ax.legend(fontsize=8)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set(xlabel="Lag (jumps)", ylabel="Bit autocorrelation", title="Bit dependence")

    fig.savefig(output / "diagnostics.png", dpi=350)
    fig.savefig(output / "diagnostics.svg")
    fig.savefig(output / "diagnostics.pdf")
    plt.close(fig)
    return output / "diagnostics.png", output / "diagnostics.pdf"


def main():
    parser = argparse.ArgumentParser(description="Plot a completed Apeiron run")
    parser.add_argument("output", type=Path, help="directory containing run results")
    args = parser.parse_args()
    for path in plot_results(args.output):
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
