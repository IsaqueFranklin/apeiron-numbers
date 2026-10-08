"""Streaming energy statistics, median digitization and three SP 800-22 tests."""

import math
import numpy as np
from scipy.special import gammaincc
from .config import EnergyConfig

BARRIER_EDGES = np.linspace(0.0, 4.0, 161)


class SampleStats:
    """Population moments merged by batches without storing individual objects."""
    def __init__(self):
        self.count = 0
        self.mean = 0.0
        self.m2 = 0.0
        self.minimum = math.inf
        self.maximum = -math.inf

    def update(self, values):
        values = np.asarray(values, dtype=float).ravel()
        if not values.size:
            return
        minimum, maximum = float(values.min()), float(values.max())
        mean = minimum if minimum == maximum else float(values.mean())
        delta = mean - self.mean
        count = self.count + values.size
        self.m2 += float(np.sum((values - mean) ** 2)) + delta**2 * self.count * values.size / count
        self.mean += delta * values.size / count
        self.count = int(count)
        self.minimum = min(self.minimum, minimum)
        self.maximum = max(self.maximum, maximum)

    def summary(self):
        return {"count": self.count, "mean_eV": self.mean,
                "std_eV": math.sqrt(self.m2 / self.count),
                "min_eV": self.minimum, "max_eV": self.maximum}


def barrier_preflight(count, seed, batch_size=100_000, config=None, samples_path=None):
    """Independent ensemble, using the exact same configuration as the landscape."""
    if count < 1 or batch_size < 1:
        raise ValueError("preflight count and batch_size must be positive")
    config = config if config is not None else EnergyConfig()
    rng = np.random.default_rng(seed)
    stats = SampleStats()
    rejected = 0
    histogram = np.zeros(len(BARRIER_EDGES) - 1, dtype=np.int64)
    saved = (np.lib.format.open_memmap(samples_path, mode="w+", dtype=np.float64, shape=(count,))
             if samples_path is not None else None)
    while stats.count < count:
        n = min(batch_size, count - stats.count)
        source = (np.full(n, config.mu_w) if config.sigma_w == 0 else
                  rng.normal(config.mu_w, config.sigma_w, n))
        target = (np.full(n, config.mu_w) if config.sigma_w == 0 else
                  rng.normal(config.mu_w, config.sigma_w, n))
        saddle = (np.full(n, config.mu_s) if config.sigma_s == 0 else
                  rng.normal(config.mu_s, config.sigma_s, n))
        bad = saddle < np.maximum(source, target)
        for iteration in range(100_000):
            if not np.any(bad):
                break
            if config.sigma_s == 0:
                raise ValueError(f"incompatible constant saddle: {config.as_dict()}")
            number = int(np.count_nonzero(bad))
            rejected += number
            saddle[bad] = rng.normal(config.mu_s, config.sigma_s, number)
            bad = saddle < np.maximum(source, target)
        else:
            raise ValueError(f"preflight rejection limit reached: {config.as_dict()}")
        barrier = saddle - source
        if not np.all(np.isfinite(barrier)):
            raise ValueError(f"nonfinite sampled energies: {config.as_dict()}")
        if saved is not None:
            saved[stats.count:stats.count+n] = barrier
        stats.update(barrier)
        histogram += np.histogram(barrier, BARRIER_EDGES)[0]
    if saved is not None:
        saved.flush()
        del saved
    return {**stats.summary(), "parameters": config.as_dict(),
            "saddle_draws": count + rejected, "saddle_rejections": rejected,
            "rejection_rate": rejected / (count + rejected),
            "untruncated_reference_mean_eV": config.mu_s - config.mu_w,
            "untruncated_reference_std_eV": math.hypot(config.sigma_w, config.sigma_s),
            "histogram": {"edges_eV": BARRIER_EDGES.tolist(), "counts": histogram.tolist(),
                          "outside_bins": count - int(histogram.sum())}}


def bits_from_waiting_times(waiting_times):
    """One bit per event, preserving order and assigning median ties to zero."""
    values = np.asarray(waiting_times, dtype=float)
    if not len(values) or not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("waiting times must be a nonempty finite nonnegative array")
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
