import math
import unittest

from apeiron.analysis import barrier_preflight, nist_subset
from apeiron.model import EnergyLandscape, FCCLattice, KMC, OFFSETS, K_B, VACANCY


class LatticeTests(unittest.TestCase):
    def test_periodic_fcc_and_composition(self):
        lattice = FCCLattice(3, 7)
        self.assertEqual(len(lattice.sites), 108)
        for site in lattice.sites:
            neighbors = lattice.neighbors(site)
            self.assertEqual(len(neighbors), 12)
            self.assertEqual(len(set(neighbors)), 12)
            self.assertTrue(all(n in lattice.occupancy for n in neighbors))
        counts = [list(lattice.occupancy.values()).count(i) for i in range(5)]
        self.assertLessEqual(max(counts) - min(counts), 1)

    def test_reversible_saddle_and_wells(self):
        lattice = FCCLattice(4, 8)
        landscape = EnergyLandscape(lattice, 9)
        origin = lattice.vacancy
        offset = OFFSETS[0]
        destination = lattice.neighbors(origin)[0]
        source_well = landscape.well(origin)
        destination_well = landscape.well(destination, (origin, destination))
        key = landscape.transition_key(origin, offset)
        saddle = landscape.saddle(origin, offset, source_well, destination_well)
        moving_species = lattice.at(destination)
        lattice.swap(destination)
        self.assertEqual(lattice.at(origin), moving_species)
        self.assertEqual(lattice.at(destination), VACANCY)
        self.assertEqual(landscape.well(destination), destination_well)
        self.assertEqual(landscape.well(origin, (destination, origin)), source_well)
        self.assertEqual(landscape.transition_key(destination, tuple(-v for v in offset)), key)
        self.assertEqual(
            landscape.saddle(destination, tuple(-v for v in offset), destination_well, source_well),
            saddle,
        )
        self.assertGreaterEqual(saddle - source_well, 0)
        self.assertGreaterEqual(saddle - destination_well, 0)
        forward = math.exp(-(saddle - source_well) / (K_B * 1273))
        backward = math.exp(-(saddle - destination_well) / (K_B * 1273))
        self.assertAlmostEqual(forward / backward, math.exp((source_well - destination_well) / (K_B * 1273)))

    def test_environment_not_position(self):
        lattice = FCCLattice(4, 10)
        landscape = EnergyLandscape(lattice, 11)
        site = lattice.vacancy
        key = landscape.well_key(site)
        value = landscape.well(site)
        self.assertEqual(value, landscape.well(site))
        neighbor = lattice.neighbors(site)[0]
        old = lattice.occupancy[neighbor]
        lattice.occupancy[neighbor] = (old + 1) % 5
        self.assertNotEqual(key, landscape.well_key(site))
        self.assertNotEqual(value, landscape.well(site))


class SimulationTests(unittest.TestCase):
    def test_event_and_reproducibility(self):
        def trajectory():
            lattice = FCCLattice(3, 12)
            landscape = EnergyLandscape(lattice, 13)
            kmc = KMC(lattice, landscape, 1273, 1e13, 14)
            return [kmc.step(i) for i in range(1, 11)]

        events = trajectory()
        self.assertEqual(events, trajectory())
        for event in events:
            self.assertEqual(len(event["barriers_eV"]), 12)
            self.assertEqual(len(event["saddles_eV"]), 12)
            self.assertEqual(len(event["rates_per_s"]), 12)
            self.assertGreater(event["dt_s"], 0)
            self.assertTrue(all(e >= 0 for e in event["barriers_eV"]))
            self.assertAlmostEqual(sum(event["rates_per_s"]), event["total_rate_per_s"])

    def test_ensemble_barriers(self):
        result = barrier_preflight(100_000, 15)
        self.assertAlmostEqual(result["mean_eV"], 0.81, delta=0.015)
        self.assertAlmostEqual(result["std_eV"], 0.32, delta=0.015)
        self.assertGreaterEqual(result["min_eV"], 0)
        self.assertGreater(result["saddle_rejections"], 0)

    def test_nist_subset(self):
        result = nist_subset([0, 1] * 100)
        self.assertEqual(result["monobit_p"], 1)
        self.assertLess(result["runs_p"], 0.01)


if __name__ == "__main__":
    unittest.main()
