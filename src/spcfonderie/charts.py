"""Cartes de contrôle : Shewhart X̄–R, règles de Western Electric, EWMA et CUSUM, capabilité.

Les fonctions de surveillance travaillent sur des moyennes réduites z = (X̄ − μ0) / (σ0 / √n)
et acceptent un tableau 1D (une série) ou 2D (plusieurs séries simulées, le temps sur le dernier axe).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Constantes des cartes X̄–R (Montgomery, Introduction to Statistical Quality Control).
CONSTANTS = {  # n: (A2, D3, D4, d2)
    2: (1.880, 0.0, 3.267, 1.128),
    3: (1.023, 0.0, 2.574, 1.693),
    4: (0.729, 0.0, 2.282, 2.059),
    5: (0.577, 0.0, 2.114, 2.326),
    6: (0.483, 0.0, 2.004, 2.534),
}


@dataclass(frozen=True)
class XbarRLimits:
    n: int
    center: float
    sigma: float          # écart-type du procédé estimé par R̄ / d2
    xbar_lcl: float
    xbar_ucl: float
    r_center: float
    r_lcl: float
    r_ucl: float

    def standardize(self, xbar: np.ndarray) -> np.ndarray:
        return (np.asarray(xbar) - self.center) / (self.sigma / np.sqrt(self.n))


def xbar_r_limits(phase1: np.ndarray) -> XbarRLimits:
    """Limites estimées sur des sous-groupes sous contrôle (phase I)."""
    n = phase1.shape[1]
    a2, d3, d4, d2 = CONSTANTS[n]
    xbarbar = float(phase1.mean())
    rbar = float(np.ptp(phase1, axis=1).mean())
    return XbarRLimits(n=n, center=xbarbar, sigma=rbar / d2,
                       xbar_lcl=xbarbar - a2 * rbar, xbar_ucl=xbarbar + a2 * rbar,
                       r_center=rbar, r_lcl=d3 * rbar, r_ucl=d4 * rbar)


def _window_sum(x: np.ndarray, w: int) -> np.ndarray:
    """Somme glissante sur les `w` derniers points (fenêtre incomplète au début)."""
    c = np.cumsum(x, axis=-1, dtype=np.int32)
    out = c.copy()
    out[..., w:] = c[..., w:] - c[..., :-w]
    return out


def western_electric(z: np.ndarray) -> dict[str, np.ndarray]:
    """Les quatre règles de Western Electric, évaluées au point qui complète le motif."""
    z = np.asarray(z)
    rules = {"1 point au-delà de 3σ": np.abs(z) > 3}
    for name, k, w, threshold in (("2 sur 3 au-delà de 2σ", 2, 3, 2.0),
                                   ("4 sur 5 au-delà de 1σ", 4, 5, 1.0),
                                   ("8 points du même côté", 8, 8, 0.0)):
        above = _window_sum((z > threshold).astype(np.int8), w) >= k
        below = _window_sum((z < -threshold).astype(np.int8), w) >= k
        rules[name] = above | below
    return rules


def ewma(z: np.ndarray, lam: float) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    s = np.zeros(z.shape[:-1])
    for t in range(z.shape[-1]):
        s = lam * z[..., t] + (1 - lam) * s
        out[..., t] = s
    return out


def ewma_limit(lam: float, L: float) -> float:
    """Limite asymptotique de l'EWMA d'une série N(0, 1)."""
    return L * np.sqrt(lam / (2 - lam))


def cusum(z: np.ndarray, k: float) -> tuple[np.ndarray, np.ndarray]:
    """CUSUM tabulaire bilatéral (C+ détecte une hausse, C− une baisse)."""
    z = np.asarray(z, dtype=float)
    hi, lo = np.empty_like(z), np.empty_like(z)
    cp = np.zeros(z.shape[:-1])
    cm = np.zeros(z.shape[:-1])
    for t in range(z.shape[-1]):
        cp = np.maximum(0.0, cp + z[..., t] - k)
        cm = np.maximum(0.0, cm - z[..., t] - k)
        hi[..., t], lo[..., t] = cp, cm
    return hi, lo


def alarms(method: str, z: np.ndarray, *, lam: float = 0.2, L: float = 3.0, k: float = 0.5,
           h: float = 5.0) -> np.ndarray:
    """Signal d'alarme de chaque méthode, point par point."""
    if method == "shewhart":
        return np.abs(z) > 3
    if method == "western_electric":
        return np.logical_or.reduce(list(western_electric(z).values()))
    if method == "ewma":
        return np.abs(ewma(z, lam)) > ewma_limit(lam, L)
    if method == "cusum":
        hi, lo = cusum(z, k)
        return (hi > h) | (lo > h)
    raise ValueError(f"méthode inconnue : {method}")


def first_alarm(flags: np.ndarray) -> np.ndarray:
    """Rang (à partir de 1) de la première alarme ; longueur de la série + 1 si aucune."""
    flags = np.asarray(flags)
    hit = flags.any(axis=-1)
    return np.where(hit, flags.argmax(axis=-1) + 1, flags.shape[-1] + 1)


def capability(x: np.ndarray, lsl: float, usl: float, sigma: float) -> tuple[float, float]:
    """Indices Cp (potentiel) et Cpk (réel, tient compte du centrage)."""
    mu = float(np.mean(x))
    cp = (usl - lsl) / (6 * sigma)
    cpk = min(usl - mu, mu - lsl) / (3 * sigma)
    return cp, cpk
