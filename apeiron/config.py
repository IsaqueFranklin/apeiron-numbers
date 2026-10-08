"""Shared, validated Gaussian landscape parameters."""

from dataclasses import asdict, dataclass
import math


@dataclass(frozen=True)
class EnergyConfig:
    mu_s: float = 0.81
    mu_w: float = 0.0
    sigma_s: float = 0.304
    sigma_w: float = 0.110

    def __post_init__(self):
        for name, value in asdict(self).items():
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite: {value}")
        if self.sigma_s < 0 or self.sigma_w < 0:
            raise ValueError(f"sigmas must be nonnegative: {asdict(self)}")
        if self.sigma_s == 0 and (self.sigma_w > 0 or self.mu_s < self.mu_w):
            raise ValueError(
                "A constant saddle cannot bound unbounded Gaussian wells, or lie "
                f"below constant wells: {asdict(self)}"
            )

    def as_dict(self):
        return {**asdict(self), "untruncated_sigma_E_eV": math.hypot(self.sigma_s, self.sigma_w)}

    def draw_saddle(self, rng, threshold):
        """Finite rejection loop; return accepted value and number of candidates."""
        if self.sigma_s == 0:
            if self.mu_s < threshold:
                raise ValueError(f"constant saddle below threshold={threshold}: {asdict(self)}")
            return self.mu_s, 1
        for draws in range(1, 100_001):
            candidate = rng.gauss(self.mu_s, self.sigma_s)
            if math.isfinite(candidate) and candidate >= threshold:
                return candidate, draws
        raise ValueError(f"saddle rejection limit reached; threshold={threshold}, parameters={asdict(self)}")
