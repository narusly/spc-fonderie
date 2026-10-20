"""Longueur moyenne de série (ARL) par simulation, et calibrage des cartes EWMA et CUSUM.

ARL0 : nombre moyen de sous-groupes avant une fausse alarme (procédé sous contrôle).
ARL1 : nombre moyen de sous-groupes avant de détecter un décalage réel.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

from .charts import alarms, first_alarm

TARGET_ARL0 = 1 / (2 * stats.norm.sf(3))   # ≈ 370,4 : la référence de la carte de Shewhart à 3σ


def shewhart_arl(shift_z: float) -> float:
    """ARL exacte de la carte de Shewhart à 3σ pour un décalage de la moyenne réduite."""
    p = stats.norm.sf(3 - shift_z) + stats.norm.cdf(-3 - shift_z)
    return 1 / p


def run_lengths(method: str, shift_z: float, n_runs: int = 4000, horizon: int = 3000,
                seed: int = 0, **params) -> np.ndarray:
    """Longueurs de série simulées. Même graine → mêmes tirages (nombres aléatoires communs)."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((n_runs, horizon)).astype(np.float32) + shift_z
    return first_alarm(alarms(method, z, **params))


def arl(method: str, shift_z: float, **kwargs) -> tuple[float, float]:
    """ARL estimée et son erreur standard Monte Carlo."""
    rl = run_lengths(method, shift_z, **kwargs)
    return float(rl.mean()), float(rl.std(ddof=1) / np.sqrt(len(rl)))


def calibrate(method: str, param: str, low: float, high: float, target: float = TARGET_ARL0,
              iterations: int = 14, **fixed) -> float:
    """Ajuste un paramètre de limite (L de l'EWMA, h du CUSUM) pour obtenir ARL0 = `target`.

    L'ARL0 croît avec la limite ; avec des nombres aléatoires communs, la recherche dichotomique converge.
    """
    for _ in range(iterations):
        mid = (low + high) / 2
        value, _ = arl(method, 0.0, **{param: mid}, **fixed)
        low, high = (mid, high) if value < target else (low, mid)
    return (low + high) / 2
