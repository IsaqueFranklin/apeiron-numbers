"""Single or multiple independently initialized scenarios with incremental output."""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gc
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import re
import shutil
import subprocess
import time

import numpy as np

from .analysis import BARRIER_EDGES, SampleStats, barrier_preflight, bits_from_waiting_times, nist_subset
from .config import EnergyConfig
from .model import EnergyLandscape, FCCLattice, KMC, SPECIES, VACANCY, K_B

FIELDS = ["event", "origin", "destination", "selected_neighbor", "moving_species", "vacancy",
          "source_well_eV", "saddles_eV", "barriers_eV", "rates_per_s", "total_rate_per_s",
          "dt_s", "time_s", "selection_u", "waiting_u"]


def write_json(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temporary.replace(path)


def provenance():
    root = Path(__file__).resolve().parents[1]
    def git(*args):
        result = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None
    files = [root / "main.py", root / "requirements.txt", *sorted((root / "apeiron").glob("*.py"))]
    return {"git_commit": git("rev-parse", "HEAD"), "git_dirty": bool(git("status", "--porcelain")),
            "source_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            "python": platform.python_version(), "platform": platform.platform(),
            "dependencies": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib")},
            "energy_key_encoding": "compact exact signatures; original repr restored for BLAKE2b seed",
            "utc_started": datetime.now(timezone.utc).isoformat()}


def parse_scenario(value):
    try:
        name, sigma_s, sigma_w = value.split(":")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name):
            raise ValueError("unsafe scenario name")
        if name.lower() in ("comparacao", "source"):
            raise ValueError("reserved scenario name")
        config = EnergyConfig(sigma_s=float(sigma_s), sigma_w=float(sigma_w))
        return name, config
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid NOME:SIGMA_S:SIGMA_W {value!r}: {exc}") from exc


def validate(args):
    for name in ("jumps", "preflight_samples", "cells"):
        if getattr(args, name) < (3 if name == "cells" else 1):
            raise ValueError(f"{name} must be at least {3 if name == 'cells' else 1}")
    for name in ("temperature", "attempt_frequency"):
        value = getattr(args, name)
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive: {value}")
    if args.progress_interval < 0 or args.seed < 0:
        raise ValueError("progress_interval and seed must be nonnegative")
    if args.block_size is not None:
        if args.block_size < 1:
            raise ValueError("block_size must be positive")
        if args.jumps >= 100 and not (
            args.block_size >= 20 and args.block_size > .01 * args.jumps
            and 1 <= args.jumps // args.block_size < 100
        ):
            raise ValueError("block_size must satisfy M >= 20, M > 0.01*n and 1 <= floor(n/M) < 100")
    if args.scenarios and (args.sigma_s is not None or args.sigma_w is not None):
        raise ValueError("use --scenarios or --sigma-s/--sigma-w, not both")
    if args.scenarios:
        scenarios = [(name, EnergyConfig(mu_s=args.mu_s, mu_w=args.mu_w,
                                        sigma_s=c.sigma_s, sigma_w=c.sigma_w)) for name, c in args.scenarios]
    else:
        scenarios = [("run", EnergyConfig(mu_s=args.mu_s, mu_w=args.mu_w,
                                         sigma_s=.304 if args.sigma_s is None else args.sigma_s,
                                         sigma_w=.110 if args.sigma_w is None else args.sigma_w))]
    if len({name for name, _ in scenarios}) != len(scenarios):
        raise ValueError("scenario names must be unique")
    return scenarios


def run_scenario(args, name, config, output, metadata):
    output.mkdir(parents=True, exist_ok=True)
    if (output / "events.csv").exists() or (output / "summary.json").exists():
        raise FileExistsError(f"refusing to overwrite existing run in {output}")
    seeds = {"composition": args.seed, "energies": args.seed+1, "kmc": args.seed+2, "preflight": args.seed+10}
    parameters = {k: v for k, v in vars(args).items() if k != "scenarios"}
    parameters.update(config.as_dict())
    parameters.update({"scenario": name, "output": str(output)})
    summary = {"status": "running", "parameters": parameters, "seeds": seeds, "provenance": metadata,
               "completed_events": 0, "requested_events": args.jumps}
    write_json(output / "summary.json", summary)
    started = time.monotonic()
    try:
        print(f"[{name}] Preflight: {args.preflight_samples:,} samples; parameters={config.as_dict()}", flush=True)
        preflight = barrier_preflight(args.preflight_samples, seeds["preflight"], config=config,
                                     samples_path=output / "preflight_barriers.npy")
        lattice = FCCLattice(args.cells, seeds["composition"])
        initial = Counter(lattice.occupancy.values())
        initial_hash = hashlib.sha256(bytes(lattice.occupancy[p] for p in lattice.sites)).hexdigest()
        landscape = EnergyLandscape(lattice, seeds["energies"], config)
        kmc = KMC(lattice, landscape, args.temperature, args.attempt_frequency, seeds["kmc"])
        dt = np.empty(args.jumps, dtype=np.float64)
        batch = np.empty((4096, 12), dtype=np.float64)
        chosen = np.empty(4096, dtype=np.float64)
        offered, selected = SampleStats(), SampleStats()
        hist_offered = np.zeros(len(BARRIER_EDGES)-1, dtype=np.int64)
        hist_selected = np.zeros_like(hist_offered)
        used = 0
        def flush_batch():
            offered.update(batch[:used])
            selected.update(chosen[:used])
            hist_offered[:] += np.histogram(batch[:used], BARRIER_EDGES)[0]
            hist_selected[:] += np.histogram(chosen[:used], BARRIER_EDGES)[0]
        kmc_started = time.monotonic()
        with (output / "events.csv").open("w", newline="", buffering=1024*1024) as events_file, \
             (output / "waiting_times.csv").open("w", newline="", buffering=1024*1024) as times_file:
            writer = csv.DictWriter(events_file, fieldnames=FIELDS)
            writer.writeheader()
            times = csv.writer(times_file)
            times.writerow(("event", "dt_s"))
            for number in range(1, args.jumps+1):
                event = kmc.step(number)
                dt[number-1] = event["dt_s"]
                batch[used] = event["barriers_eV"]
                chosen[used] = event["barriers_eV"][event["selected_neighbor"]]
                used += 1
                writer.writerow({key: json.dumps(value, separators=(",", ":")) if isinstance(value, (list, tuple))
                                 else value for key, value in event.items()})
                times.writerow((number, event["dt_s"]))
                summary["completed_events"] = number
                if used == len(batch):
                    flush_batch()
                    used = 0
                if args.progress_interval and (number % args.progress_interval == 0 or number == args.jumps):
                    elapsed = time.monotonic()-kmc_started
                    rate = number/elapsed
                    print(f"[{name}] {number:,}/{args.jumps:,} ({number/args.jumps:.1%}) | "
                          f"elapsed {elapsed:.0f}s | {rate:.0f} jumps/s | ETA {(args.jumps-number)/rate:.0f}s", flush=True)
                    events_file.flush()
                    times_file.flush()
                    write_json(output / "summary.json", summary)
            if used:
                flush_batch()
        kmc_duration = time.monotonic()-kmc_started
        final = Counter(lattice.occupancy.values())
        if initial != final or final[VACANCY] != 1:
            raise RuntimeError("species/vacancy conservation failed")
        bits, threshold = bits_from_waiting_times(dt)
        (output / "bits.txt").write_bytes((bits + 48).tobytes() + b"\n")
        np.savez(output / "barrier_histograms.npz", edges_eV=BARRIER_EDGES,
                 offered_counts=hist_offered, selected_counts=hist_selected)
        probs = [.01, .05, .25, .75, .95, .99]
        waiting = {"count": len(dt), "mean_s": float(dt.mean()), "median_s": threshold,
                   "std_s": float(dt.std()), "min_s": float(dt.min()), "max_s": float(dt.max()),
                   "zero_count": int(np.count_nonzero(dt == 0)),
                   "quantiles_s": {f"{int(p*100):02d}%": float(q) for p, q in zip(probs, np.quantile(dt, probs))}}
        summary.update({"preflight": preflight, "realized_offered_barriers": offered.summary(),
                        "selected_barriers": selected.summary(),
                        "landscape": {"unique_wells": len(landscape.wells), "unique_saddles": len(landscape.saddles),
                                      "saddle_draws": landscape.saddle_draws, "saddle_rejections": landscape.saddle_rejections,
                                      "rejection_rate": landscape.saddle_rejections/landscape.saddle_draws},
                        "lattice": {"sites": len(lattice.sites), "atoms": len(lattice.sites)-1, "vacancies": final[VACANCY],
                                    "initial_composition": {s: initial[i] for i, s in enumerate(SPECIES)},
                                    "final_composition": {s: final[i] for i, s in enumerate(SPECIES)},
                                    "initial_occupancy_sha256": initial_hash},
                        "waiting_times": waiting, "total_time_s": kmc.time,
                        "digitization": {"rule": "dt > whole-run median", "median_dt_s": threshold,
                                         "ties": int(np.count_nonzero(dt == threshold)), "bits": len(bits),
                                         "ones_fraction": float(bits.mean())},
                        "nist_sp800_22_subset": nist_subset(bits, block_size=args.block_size), "kmc_duration_s": kmc_duration})
        if config.sigma_s == config.sigma_w == 0:
            r0 = 12*args.attempt_frequency*math.exp(-(config.mu_s-config.mu_w)/(K_B*args.temperature))
            summary["homogeneous_reference"] = {"total_rate_per_s": r0, "mean_waiting_s": 1/r0,
                                                "observed_mean_relative_error": waiting["mean_s"]*r0-1}
        # Release large caches before loading plot data or initializing another run.
        del kmc, landscape, lattice
        gc.collect()
        if args.plot:
            from .plotting import plot_results
            plot_results(output, dt=dt, bits=bits)
        summary.update({"status": "complete", "duration_s": time.monotonic()-started,
                        "utc_finished": datetime.now(timezone.utc).isoformat()})
        write_json(output / "summary.json", summary)
        print(f"[{name}] Complete: {args.jumps:,} events/bits in {summary['duration_s']:.1f}s → {output}", flush=True)
        return summary
    except BaseException as exc:
        summary.update({"status": "incomplete", "duration_s": time.monotonic()-started,
                        "failure": f"{type(exc).__name__}: {exc}"})
        write_json(output / "summary.json", summary)
        raise


def run(args):
    scenarios = validate(args)
    root = Path(args.output)
    root.mkdir(parents=True, exist_ok=True)
    for name, _ in scenarios:
        directory = root / name if args.scenarios else root
        if (directory / "events.csv").exists() or (directory / "summary.json").exists():
            raise FileExistsError(f"refusing to overwrite existing run in {directory}")
    estimated = sum(args.jumps*(1000 if config.sigma_s == config.sigma_w == 0 else 1600)
                    + args.preflight_samples*8 for _, config in scenarios) + 30_000_000
    free = shutil.disk_usage(root).free
    if free < estimated:
        raise OSError(f"insufficient storage: {free/1e9:.2f} GB free; conservative estimate {estimated/1e9:.2f} GB "
                      "for complete event records. Free disk space or choose another --output volume.")
    metadata = provenance()
    multi = bool(args.scenarios)
    if multi and (root / "manifest.json").exists():
        raise FileExistsError(f"manifest already exists: {root}")
    source = root / "source"
    source.mkdir(exist_ok=True)
    project = Path(__file__).resolve().parents[1]
    for relative in metadata["source_sha256"]:
        target = source / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(project / relative, target)
    manifest = {"status": "running", "provenance": metadata, "argv": list(__import__("sys").argv),
                "storage": {"initial_free_bytes": free, "estimated_bytes": estimated},
                "seeds": {"composition": args.seed, "energies": args.seed+1, "kmc": args.seed+2, "preflight": args.seed+10},
                "scenarios": [{"name": name, "parameters": config.as_dict(), "status": "pending"} for name, config in scenarios]}
    write_json(root / "manifest.json", manifest)
    started = time.monotonic()
    try:
        for entry, (name, config) in zip(manifest["scenarios"], scenarios):
            entry["status"] = "running"
            write_json(root / "manifest.json", manifest)
            result = run_scenario(args, name, config, root / name if multi else root, metadata)
            entry.update({"status": "complete", "duration_s": result["duration_s"],
                          "summary": f"{name}/summary.json" if multi else "summary.json"})
            write_json(root / "manifest.json", manifest)
        if multi:
            from .comparison import compare_results
            compare_results(root, plot=args.plot)
        manifest.update({"status": "complete", "duration_s": time.monotonic()-started,
                         "utc_finished": datetime.now(timezone.utc).isoformat()})
        write_json(root / "manifest.json", manifest)
    except BaseException as exc:
        manifest.update({"status": "incomplete", "failure": f"{type(exc).__name__}: {exc}",
                         "duration_s": time.monotonic()-started})
        for entry in manifest["scenarios"]:
            if entry["status"] == "running":
                entry["status"] = "incomplete"
        write_json(root / "manifest.json", manifest)
        raise
    return manifest


def build_parser():
    parser = argparse.ArgumentParser(description="Vacancy diffusion on a periodic FCC alloy")
    parser.add_argument("--cells", type=int, default=8, help="conventional FCC cells per axis")
    parser.add_argument("--jumps", type=int, default=10_000)
    parser.add_argument("--temperature", type=float, default=1273.0, help="kelvin")
    parser.add_argument("--attempt-frequency", type=float, default=1e13, help="per second")
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--block-size", type=int, default=None,
                        help="SP 800-22 block size; default max(20, ceil(jumps/50)), 20000 for 1M bits")
    parser.add_argument("--preflight-samples", type=int, default=1_000_000)
    parser.add_argument("--mu-s", type=float, default=.81, help="mean saddle energy, eV")
    parser.add_argument("--mu-w", type=float, default=0.0, help="mean well energy, eV")
    parser.add_argument("--sigma-s", type=float, default=None, help="single-run saddle sigma, eV (default .304)")
    parser.add_argument("--sigma-w", type=float, default=None, help="single-run well sigma, eV (default .110)")
    parser.add_argument("--scenarios", type=parse_scenario, nargs="+", metavar="NAME:SIGMA_S:SIGMA_W")
    parser.add_argument("--progress-interval", type=int, default=10_000, help="progress every N jumps; 0 disables it")
    parser.add_argument("--output", default="results/run")
    parser.add_argument("--plot", action="store_true")
    return parser


def main():
    parser = build_parser()
    try:
        run(parser.parse_args())
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
