"""Barrier calibration, exploratory bit extraction and statistical summaries."""

import math

import numpy as np
from scipy.special import gammaincc


def barrier_preflight(count, seed, batch_size=100_000):
    """Draw independent model environments and reject saddles below both wells.

    This is an ensemble check, not the distribution visited by a KMC walk.
    """
    if count < 1:
        raise ValueError("preflight count must be positive")
    rng = np.random.default_rng(seed)
    accepted = 0
    rejected = 0
    total = 0.0
    total_sq = 0.0
    minimum = math.inf
    while accepted < count:
        n = min(batch_size, count - accepted)
        source = rng.normal(0, 0.110, n)
        target = rng.normal(0, 0.110, n)
        saddle = rng.normal(0.81, 0.304, n)
        bad = saddle < np.maximum(source, target)
        while np.any(bad):
            rejected += int(np.count_nonzero(bad))
            saddle[bad] = rng.normal(0.81, 0.304, int(np.count_nonzero(bad)))
            bad = saddle < np.maximum(source, target)
        barrier = saddle - source
        total += float(np.sum(barrier))
        total_sq += float(np.dot(barrier, barrier))
        minimum = min(minimum, float(np.min(barrier)))
        accepted += n
    mean = total / count
    return {
        "count": count, "mean_eV": mean,
        "std_eV": math.sqrt(max(0, total_sq / count - mean * mean)),
        "min_eV": minimum, "saddle_draws": count + rejected,
        "saddle_rejections": rejected,
        "rejection_rate": rejected / (count + rejected),
        "untruncated_reference_mean_eV": 0.81,
        "untruncated_reference_std_eV": math.hypot(0.110, 0.304),
    }


def bits_from_waiting_times(waiting_times):
    """Exploratory median digitization; the median is estimated on this run."""
    values = np.asarray(waiting_times, dtype=float)
    if len(values) == 0:
        raise ValueError("no waiting times")
    threshold = float(np.median(values))
    return (values > threshold).astype(np.uint8), threshold


def nist_subset(bits, block_size=128):
    """SP 800-22 Rev. 1a monobit, block frequency, and runs tests."""
    x = np.asarray(bits, dtype=np.uint8)
    n = len(x)
    if n < 100:
        return {"status": "insufficient_bits", "bits": n, "minimum_bits": 100}
    ones = int(x.sum())
    proportion = ones / n
    monobit = math.erfc(abs(2 * ones - n) / math.sqrt(2 * n))
    blocks = n // block_size
    if blocks:
        chunk = x[: blocks * block_size].reshape(blocks, block_size)
        chi_square = float(4 * block_size * np.sum((chunk.mean(axis=1) - 0.5) ** 2))
        block_p = float(gammaincc(blocks / 2, chi_square / 2))
    else:
        block_p = None
    if abs(proportion - 0.5) >= 2 / math.sqrt(n):
        runs_p = None
    else:
        runs = 1 + int(np.count_nonzero(x[1:] != x[:-1]))
        numerator = abs(runs - 2 * n * proportion * (1 - proportion))
        denominator = 2 * math.sqrt(2 * n) * proportion * (1 - proportion)
        runs_p = math.erfc(numerator / denominator)
    alpha = 0.01
    return {
        "status": "ok", "bits": n, "ones_fraction": proportion,
        "alpha": alpha, "block_size": block_size,
        "monobit_p": monobit, "block_frequency_p": block_p,
        "runs_p": runs_p,
        "runs_prerequisite_met": runs_p is not None,
        "passed": {
            "monobit": monobit >= alpha,
            "block_frequency": block_p >= alpha if block_p is not None else None,
            "runs": runs_p >= alpha if runs_p is not None else None,
        },
    }
