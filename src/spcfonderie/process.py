"""Simulation du sable de fonderie : trois caractéristiques mesurées en sous-groupes, avec défauts injectés.

Les cibles et dispersions reprennent des ordres de grandeur courants en fonderie (sable à vert) ;
les tolérances sont hypothétiques et servent au calcul de capabilité.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Characteristic:
    label: str
    unit: str
    target: float
    sigma: float
    lsl: float
    usl: float


CHARACTERISTICS = {
    "humidite": Characteristic("Humidité", "%", 4.5, 0.28, 3.6, 5.4),
    "granulometrie": Characteristic("Granulométrie", "indice AFS", 55.0, 3.8, 43.0, 67.0),
    "resistance": Characteristic("Résistance à la compression", "kPa", 120.0, 9.7, 90.0, 150.0),
}


@dataclass(frozen=True)
class Fault:
    """Défaut à partir du sous-groupe `start`.

    - shift : la moyenne se décale de `size` écarts-types ;
    - drift : la moyenne dérive de `size` écart-type par sous-groupe ;
    - variance : l'écart-type est multiplié par `size`.
    """
    kind: str
    start: int
    size: float


def simulate(c: Characteristic, n_subgroups: int, n: int, rng: np.random.Generator,
             fault: Fault | None = None) -> np.ndarray:
    """Mesures individuelles, tableau (sous-groupes × n)."""
    mean = np.full(n_subgroups, c.target)
    sigma = np.full(n_subgroups, c.sigma)
    if fault is not None:
        after = np.arange(n_subgroups) >= fault.start
        if fault.kind == "shift":
            mean[after] += fault.size * c.sigma
        elif fault.kind == "drift":
            mean[after] += fault.size * c.sigma * (np.arange(n_subgroups)[after] - fault.start + 1)
        elif fault.kind == "variance":
            sigma[after] *= fault.size
        else:
            raise ValueError(f"défaut inconnu : {fault.kind}")
    return mean[:, None] + sigma[:, None] * rng.standard_normal((n_subgroups, n))
