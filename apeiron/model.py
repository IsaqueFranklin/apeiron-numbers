"""Periodic FCC lattice and rejection-free residence-time KMC."""

import hashlib
import math
import random

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

    def at(self, point):
        return self.occupancy[tuple(v % self.side for v in point)]

    def neighbors(self, point):
        return tuple(tuple((a + b) % self.side for a, b in zip(point, d)) for d in OFFSETS)

    def swap(self, destination):
        origin = self.vacancy
        if destination not in self.neighbors(origin):
            raise ValueError("destination is not a nearest neighbor")
        self.occupancy[origin] = self.occupancy[destination]
        self.occupancy[destination] = VACANCY
        self.vacancy = destination


class EnergyLandscape:
    """Translation-invariant chemical keys with reproducible Gaussian energies."""

    def __init__(self, lattice, seed):
        self.lattice = lattice
        self.seed = seed
        self.wells = {}
        self.saddles = {}
        self.saddle_draws = 0
        self.saddle_rejections = 0

    def _rng(self, kind, key):
        material = repr((self.seed, kind, key)).encode("ascii")
        digest = hashlib.blake2b(material, digest_size=16).digest()
        return random.Random(int.from_bytes(digest, "big"))

    def well_key(self, site, virtual_swap=None):
        """Ordered 12-neighbor chemical shell; no absolute coordinates in key."""
        def atom(point):
            if virtual_swap is not None:
                a, b = virtual_swap
                if point == a:
                    return self.lattice.at(b)
                if point == b:
                    return self.lattice.at(a)
            return self.lattice.at(point)

        return tuple(atom(p) for p in self.lattice.neighbors(site))

    def well(self, site, virtual_swap=None):
        key = self.well_key(site, virtual_swap)
        if key not in self.wells:
            self.wells[key] = self._rng("well", key).gauss(0.0, 0.110)
        return self.wells[key]

    def transition_key(self, origin, offset):
        """Moving atom plus bridge shell, invariant to reversing the hop."""
        destination = tuple((a + b) % self.lattice.side for a, b in zip(origin, offset))
        moving_atom = self.lattice.at(destination)
        bridge = set(OFFSETS)
        bridge.update(tuple(a + b for a, b in zip(offset, d)) for d in OFFSETS)
        bridge.update(((0, 0, 0), offset))

        def orientation(anchor, coordinates, edge_offset):
            return tuple(
                (relative, self.lattice.at(
                    tuple((a + b) % self.lattice.side for a, b in zip(anchor, relative))
                ))
                for relative in sorted(coordinates)
                if relative != (0, 0, 0) and relative != edge_offset
            )

        # Endpoint atoms are omitted from the bridge shell. The atom moving
        # through the saddle is recorded separately. Reverse orientation maps
        # every relative vector v to v - offset.
        forward = orientation(origin, bridge, offset)
        opposite = tuple(-v for v in offset)
        reverse = orientation(destination, {tuple(v - d for v, d in zip(p, offset)) for p in bridge}, opposite)
        return (moving_atom, min(forward, reverse))

    def saddle(self, origin, offset, source_well, destination_well):
        key = self.transition_key(origin, offset)
        if key not in self.saddles:
            rng = self._rng("saddle", key)
            threshold = max(source_well, destination_well)
            while True:
                candidate = rng.gauss(0.81, 0.304)
                self.saddle_draws += 1
                if candidate >= threshold:
                    self.saddles[key] = candidate
                    break
                self.saddle_rejections += 1
        saddle = self.saddles[key]
        if saddle < max(source_well, destination_well):
            raise RuntimeError("same transition key has inconsistent endpoint wells")
        return saddle


class KMC:
    def __init__(self, lattice, landscape, temperature, attempt_frequency, seed):
        if temperature <= 0 or attempt_frequency <= 0:
            raise ValueError("temperature and attempt frequency must be positive")
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
