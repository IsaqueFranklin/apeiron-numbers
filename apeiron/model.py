"""Periodic FCC lattice and rejection-free residence-time KMC."""

import hashlib
import math
import random

from .config import EnergyConfig

K_B = 8.617333262e-5  # eV/K
SPECIES = ("Co", "Ni", "Cr", "Fe", "Mn")
VACANCY = 5
# FCC represented by integer coordinates with even x+y+z parity.
OFFSETS = tuple(
    (x, y, z)
    for x in (-1, 0, 1)
    for y in (-1, 0, 1)
    for z in (-1, 0, 1)
    if x * x + y * y + z * z == 2
)


class FCCLattice:
    """A conventional-cell periodic FCC box with exactly one vacancy."""

    def __init__(self, cells, seed):
        if cells < 3:
            raise ValueError("cells must be at least 3 to avoid duplicate neighbors")
        self.side = 2 * cells
        self.sites = tuple(
            (x, y, z)
            for x in range(self.side)
            for y in range(self.side)
            for z in range(self.side)
            if (x + y + z) % 2 == 0
        )
        rng = random.Random(seed)
        self.vacancy = self.sites[0]
        atoms = [i % 5 for i in range(len(self.sites) - 1)]
        rng.shuffle(atoms)
        self.occupancy = dict(zip(self.sites[1:], atoms))
        self.occupancy[self.vacancy] = VACANCY
        self._neighbors = {
            site: tuple(tuple((a + b) % self.side for a, b in zip(site, d)) for d in OFFSETS)
            for site in self.sites
        }

    def at(self, point):
        return self.occupancy[tuple(v % self.side for v in point)]

    def neighbors(self, point):
        point = tuple(point)
        try:
            return self._neighbors[point]
        except KeyError:
            return self._neighbors[tuple(v % self.side for v in point)]

    def swap(self, destination):
        origin = self.vacancy
        if destination not in self.neighbors(origin):
            raise ValueError("destination is not a nearest neighbor")
        self.occupancy[origin] = self.occupancy[destination]
        self.occupancy[destination] = VACANCY
        self.vacancy = destination


class EnergyLandscape:
    """Translation-invariant chemical keys with reproducible Gaussian energies."""

    def __init__(self, lattice, seed, config=None):
        self.lattice = lattice
        self.seed = seed
        self.config = config if config is not None else EnergyConfig()
        self.wells = {}
        self.saddles = {}
        self.saddle_draws = 0
        self.saddle_rejections = 0
        self._frames = {}
        self._bridge_sites = {}
        for offset in OFFSETS:
            bridge = set(OFFSETS)
            bridge.update(tuple(a + b for a, b in zip(offset, d)) for d in OFFSETS)
            forward = tuple(sorted(bridge - {(0, 0, 0), offset}))
            reverse = tuple(sorted(tuple(v - d for v, d in zip(p, offset)) for p in forward))
            backwards = reverse < forward
            canonical = tuple(-v for v in offset) if backwards else offset
            self._frames[offset] = (OFFSETS.index(canonical), min(forward, reverse), backwards)

    def _rng(self, kind, key):
        # Decode compact, exact cache keys to preserve the original seeded draws.
        if kind == "well":
            original_key = tuple(key)
        else:
            coordinates = self._frames[OFFSETS[key[0]]][1]
            original_key = (key[1], tuple(zip(coordinates, key[2:])))
        material = repr((self.seed, kind, original_key)).encode("ascii")
        digest = hashlib.blake2b(material, digest_size=16).digest()
        return random.Random(int.from_bytes(digest, "big"))

    def well_key(self, site, virtual_swap=None):
        """Ordered 12-neighbor chemical shell; no absolute coordinates in key."""
        def atom(point):
            if virtual_swap is not None:
                a, b = virtual_swap
                if point == a:
                    return self.lattice.occupancy[b]
                if point == b:
                    return self.lattice.occupancy[a]
            return self.lattice.occupancy[point]

        return bytes(atom(p) for p in self.lattice.neighbors(site))

    def well(self, site, virtual_swap=None):
        key = self.well_key(site, virtual_swap)
        if key not in self.wells:
            self.wells[key] = (self.config.mu_w if self.config.sigma_w == 0 else
                               self._rng("well", key).gauss(self.config.mu_w, self.config.sigma_w))
        return self.wells[key]

    def transition_key(self, origin, offset):
        """Moving atom plus bridge shell, invariant to reversing the hop."""
        destination = self.lattice.neighbors(origin)[OFFSETS.index(offset)]
        direction, coordinates, backwards = self._frames[offset]
        anchor = destination if backwards else origin
        geometry_key = (anchor, direction)
        if geometry_key not in self._bridge_sites:
            self._bridge_sites[geometry_key] = tuple(
                tuple((a + b) % self.lattice.side for a, b in zip(anchor, p))
                for p in coordinates
            )
        occupancy = self.lattice.occupancy
        return bytes((direction, occupancy[destination])) + bytes(
            occupancy[p] for p in self._bridge_sites[geometry_key]
        )

    def saddle(self, origin, offset, source_well, destination_well):
        key = self.transition_key(origin, offset)
        if key not in self.saddles:
            threshold = max(source_well, destination_well)
            rng = self._rng("saddle", key) if self.config.sigma_s else None
            candidate, draws = self.config.draw_saddle(rng, threshold)
            self.saddles[key] = candidate
            self.saddle_draws += draws
            self.saddle_rejections += draws - 1
        saddle = self.saddles[key]
        if saddle < max(source_well, destination_well):
            raise RuntimeError("same transition key has inconsistent endpoint wells")
        return saddle


class KMC:
    def __init__(self, lattice, landscape, temperature, attempt_frequency, seed):
        if not all(math.isfinite(v) and v > 0 for v in (temperature, attempt_frequency)):
            raise ValueError("temperature and attempt frequency must be finite and positive")
        self.lattice = lattice
        self.landscape = landscape
        self.temperature = temperature
        self.attempt_frequency = attempt_frequency
        self.rng = random.Random(seed)
        self.time = 0.0

    def step(self, number):
        origin = self.lattice.vacancy
        source_well = self.landscape.well(origin)
        destinations = self.lattice.neighbors(origin)
        saddles = []
        barriers = []
        rates = []
        for offset, destination in zip(OFFSETS, destinations):
            target_well = self.landscape.well(destination, (origin, destination))
            saddle = self.landscape.saddle(origin, offset, source_well, target_well)
            barrier = saddle - source_well
            saddles.append(saddle)
            barriers.append(barrier)
            rates.append(self.attempt_frequency * math.exp(-barrier / (K_B * self.temperature)))
        total_rate = math.fsum(rates)
        if not math.isfinite(total_rate) or total_rate <= 0:
            raise ArithmeticError("total rate is zero or not finite; adjust parameters")
        selection_u = self.rng.random()
        waiting_u = 1.0 - self.rng.random()  # U(0, 1], including 1
        target = selection_u * total_rate
        cumulative = 0.0
        selected = 11
        for index, rate in enumerate(rates):
            cumulative += rate
            if target < cumulative:
                selected = index
                break
        dt = -math.log(waiting_u) / total_rate
        self.time += dt
        destination = destinations[selected]
        moving_species = SPECIES[self.lattice.at(destination)]
        record = {
            "event": number, "origin": origin, "destination": destination,
            "selected_neighbor": selected, "moving_species": moving_species,
            "vacancy": destination, "source_well_eV": source_well,
            "saddles_eV": saddles, "barriers_eV": barriers,
            "rates_per_s": rates, "total_rate_per_s": total_rate,
            "dt_s": dt, "time_s": self.time,
            "selection_u": selection_u, "waiting_u": waiting_u,
        }
        self.lattice.swap(destination)
        return record
