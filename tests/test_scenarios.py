import csv
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from collections import Counter

import numpy as np

from apeiron.analysis import barrier_preflight, bits_from_waiting_times, nist_subset, select_block_size
from apeiron.cli import build_parser, run, validate
from apeiron.config import EnergyConfig
from apeiron.model import FCCLattice, EnergyLandscape, KMC, K_B, OFFSETS


class ScientificScenarioTests(unittest.TestCase):
    def test_parameter_validation(self):
        for params in ({"sigma_s": -1}, {"sigma_w": float("nan")}, {"mu_s": float("inf")},
                       {"sigma_s": 0, "sigma_w": .1}, {"sigma_s": 0, "sigma_w": 0, "mu_s": -.1}):
            with self.subTest(params=params), self.assertRaises(ValueError):
                EnergyConfig(**params)
        for flag, value in (("--temperature", "nan"), ("--attempt-frequency", "inf"),
                            ("--cells", "2"), ("--jumps", "0"), ("--preflight-samples", "0")):
            with self.subTest(flag=flag), self.assertRaises(ValueError):
                validate(build_parser().parse_args([flag, value]))

    def test_parameter_propagation_and_homogeneous_limit(self):
        config = EnergyConfig(sigma_s=0, sigma_w=0)
        landscape = EnergyLandscape(FCCLattice(8, 2020), 2021, config)
        self.assertIs(landscape.config, config)
        check = barrier_preflight(1000, 2030, config=config)
        self.assertEqual(check["parameters"], config.as_dict())
        self.assertEqual(check["mean_eV"], .81)
        self.assertEqual(check["std_eV"], 0)
        self.assertEqual(check["saddle_rejections"], 0)
        kmc = KMC(landscape.lattice, landscape, 1273, 1e13, 2022)
        expected = 12*1e13*math.exp(-.81/(K_B*1273))
        event = kmc.step(1)
        self.assertEqual(event["barriers_eV"], [.81]*12)
        self.assertAlmostEqual(event["total_rate_per_s"]/expected, 1)
        self.assertEqual(event["dt_s"], -math.log(event["waiting_u"])/event["total_rate_per_s"])
        self.assertAlmostEqual(1/expected, 1.34150e-11, delta=1e-16)
        h1 = barrier_preflight(100000, 2030, config=EnergyConfig(sigma_s=.152, sigma_w=.055))
        self.assertAlmostEqual(h1["std_eV"], math.hypot(.152,.055), delta=.002)
        self.assertNotEqual(h1["parameters"], check["parameters"])

    def test_compact_keys_preserve_original_energies_and_trajectory(self):
        # Independent implementation at the scientific reference commit.
        source = subprocess.check_output(["git", "show", "f50cda4fd5a7c8446a8676460a3d50550b29919a:apeiron/model.py"], text=True)
        historical = {}
        exec(compile(source, "historical-model", "exec"), historical)
        old_lattice = historical["FCCLattice"](3, 2020)
        old = historical["KMC"](old_lattice, historical["EnergyLandscape"](old_lattice,2021),1273,1e13,2022)
        lattice = FCCLattice(3,2020)
        new = KMC(lattice, EnergyLandscape(lattice,2021),1273,1e13,2022)
        for i in range(1,101):
            self.assertEqual(new.step(i), old.step(i))

    def test_conservation_and_separate_caches(self):
        lattice = FCCLattice(8,2020)
        before = Counter(lattice.occupancy.values())
        self.assertEqual([before[i] for i in range(6)], [410,410,409,409,409,1])
        landscapes = [EnergyLandscape(FCCLattice(3,2020),2021,EnergyConfig(sigma_s=s,sigma_w=w))
                      for s,w in ((0,0),(.152,.055),(.304,.110))]
        self.assertIsNot(landscapes[0].wells, landscapes[1].wells)
        self.assertIsNot(landscapes[1].saddles, landscapes[2].saddles)
        kmc = KMC(lattice, EnergyLandscape(lattice,2021),1273,1e13,2022)
        for i in range(300): kmc.step(i)
        self.assertEqual(Counter(lattice.occupancy.values()), before)

    def test_translation_and_reversal_all_directions(self):
        for offset in OFFSETS:
            lattice = FCCLattice(3,2020)
            self.assertEqual(lattice.neighbors((lattice.side,0,0)), lattice.neighbors((0,0,0)))
            landscape = EnergyLandscape(lattice,2021)
            origin = lattice.vacancy
            dest = lattice.neighbors(origin)[OFFSETS.index(offset)]
            a = landscape.well(origin)
            b = landscape.well(dest,(origin,dest))
            key = landscape.transition_key(origin,offset)
            saddle = landscape.saddle(origin,offset,a,b)
            lattice.swap(dest)
            self.assertEqual(landscape.transition_key(dest,tuple(-v for v in offset)),key)
            self.assertEqual(landscape.saddle(dest,tuple(-v for v in offset),b,a),saddle)
            shifted = FCCLattice(3,2020)
            shift = (2,0,0)
            transform = lambda p: tuple((a+b)%lattice.side for a,b in zip(p,shift))
            shifted.occupancy = {transform(p): value for p,value in lattice.occupancy.items()}
            other = EnergyLandscape(shifted,2021)
            self.assertEqual(other.well(transform(dest)),landscape.well(dest))
            self.assertEqual(other.transition_key(transform(dest),tuple(-v for v in offset)),key)

    def test_ties_and_not_applicable(self):
        bits, threshold = bits_from_waiting_times([1,2,2,3])
        self.assertEqual(threshold,2)
        self.assertEqual(bits.tolist(),[0,0,0,1])
        self.assertEqual(nist_subset([0]*200)["tests"]["runs"]["status"],"not_applicable")
        short = nist_subset([0,1]*20)
        self.assertIsNone(short["monobit_p"])
        self.assertEqual(short["tests"]["monobit"]["status"],"not_applicable")
        no_blocks = nist_subset([0,1]*50,block_size=128)
        self.assertIsNone(no_blocks["block_frequency_p"])
        million = nist_subset(np.tile([0,1],500000))
        self.assertEqual((million["block_size"],million["complete_blocks"],million["residual_bits"]),(20000,50,0))
        self.assertTrue(all(million["block_recommendations"].values()))
        self.assertEqual(million["runs_p"],0)
        self.assertEqual(million["tests"]["runs"]["status"],"fail")

    def test_nist_block_recommendations_and_reference_example(self):
        for n in (100, 199, 1000, 20001, 999999, 1000000):
            m = select_block_size(n)
            self.assertGreaterEqual(m, 20)
            self.assertGreater(m, .01*n)
            self.assertGreaterEqual(n//m, 1)
            self.assertLess(n//m, 100)
        # Published SP 800-22 section 2.2.8 example, independent of auto sizing.
        text = ("11001001000011111101101010100010001000010110100011"
                "00001000110100110001001100011001100010100010111000")
        reference = nist_subset([int(c) for c in text], block_size=10)
        self.assertAlmostEqual(reference['tests']['block_frequency']['statistic'], 7.2)
        self.assertAlmostEqual(reference['block_frequency_p'], .706438, places=6)
        clustered = nist_subset(np.r_[np.zeros(500000, dtype=np.uint8), np.ones(500000, dtype=np.uint8)])
        self.assertEqual(clustered['monobit_p'], 1)
        self.assertEqual(clustered['tests']['block_frequency']['statistic'], 1000000)
        self.assertEqual(clustered['block_frequency_p'], 0)
        self.assertEqual(clustered['runs_count'], 2)
        for value in (0, -1, 1.5, True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                nist_subset([0,1]*100, block_size=value)

    def test_cli_block_size_validation(self):
        for m in (0, -1, 128, 10000):
            with self.subTest(m=m), self.assertRaises(ValueError):
                validate(build_parser().parse_args(['--jumps','1000000','--block-size',str(m)]))
        validate(build_parser().parse_args(['--jumps','1000000','--block-size','20000']))

    def test_cli_three_scenarios_correspondence_and_reproducibility(self):
        with tempfile.TemporaryDirectory(prefix="apeiron-test-") as temp:
            all_summaries = []
            for repetition in ("a","b"):
                output = Path(temp)/repetition
                args = build_parser().parse_args(["--scenarios","H0:0:0","H1:0.152:0.055","H2:0.304:0.110",
                    "--jumps","200","--preflight-samples","1000","--cells","3","--progress-interval","0","--output",str(output)])
                manifest = run(args)
                self.assertEqual(manifest["status"],"complete")
                preserved = (output / "manifest.json").read_bytes()
                with self.assertRaises(FileExistsError):
                    run(args)
                self.assertEqual((output / "manifest.json").read_bytes(), preserved)
                summaries = []
                for name in ("H0","H1","H2"):
                    folder = output/name
                    summary = json.loads((folder/"summary.json").read_text())
                    self.assertEqual(summary['nist_sp800_22_subset']['block_size'], 20)
                    self.assertTrue(all(summary['nist_sp800_22_subset']['block_recommendations'].values()))
                    with (folder/"events.csv").open() as handle: events=list(csv.DictReader(handle))
                    with (folder/"waiting_times.csv").open() as handle: waiting=list(csv.DictReader(handle))
                    bits = (folder/"bits.txt").read_text().strip()
                    self.assertEqual((len(events),len(waiting),len(bits)),(200,200,200))
                    self.assertEqual(summary["realized_offered_barriers"]["count"],2400)
                    self.assertEqual(summary["selected_barriers"]["count"],200)
                    for event,row,bit in zip(events,waiting,bits):
                        self.assertEqual(event["event"],row["event"])
                        self.assertEqual(event["dt_s"],row["dt_s"])
                        self.assertEqual(int(bit),int(float(row["dt_s"]) > summary["digitization"]["median_dt_s"]))
                    summaries.append(summary)
                self.assertEqual(len({s["lattice"]["initial_occupancy_sha256"] for s in summaries}),1)
                all_summaries.append(summaries)
            for name,a,b in zip(("H0","H1","H2"),*all_summaries):
                self.assertEqual(a["waiting_times"],b["waiting_times"])
                self.assertEqual((Path(temp)/"a"/name/"events.csv").read_bytes(),(Path(temp)/"b"/name/"events.csv").read_bytes())


if __name__ == "__main__": unittest.main()
