"""Analyse complète : `python -m spcfonderie.report` écrit `reports/resultats.md` et les figures."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import style
from .arl import TARGET_ARL0, arl, calibrate, run_lengths, shewhart_arl
from .charts import alarms, capability, first_alarm, xbar_r_limits
from .process import CHARACTERISTICS, Fault, simulate

ROOT = Path(__file__).resolve().parents[2]
REPORTS = ROOT / "reports"
FIGURES = REPORTS / "figures"
N, PHASE1, PHASE2, SEED = 4, 30, 100, 7
LAMBDA, K = 0.2, 0.5
# Scénario hypothétique de phase II : un défaut différent sur chaque caractéristique.
SCENARIO = {
    "humidite": Fault("shift", 40, 0.5),        # dosage d'eau déréglé : +0,5σ d'un coup
    "resistance": Fault("drift", 50, -0.05),    # dérive lente : −0,05σ par sous-groupe
    "granulometrie": Fault("variance", 60, 2.0),  # dispersion doublée
}


METHODS = {"shewhart": "Shewhart 3σ", "western_electric": "Shewhart + règles WE",
           "ewma": "EWMA (λ = 0,2)", "cusum": "CUSUM (k = 0,5)"}
COLORS = dict(zip(METHODS, style.SERIES))      # une couleur par méthode, la même dans toutes les figures
DETECTION_SHIFT = 0.5                          # décalage étudié en distribution (en σ du procédé)


def signed(x: float, decimals: int) -> str:
    return f"{x:+.{decimals}f}".replace(".", ",").replace("-", "−")


def fault_label(f: Fault) -> str:
    if f.kind == "shift":
        return f"décalage de {signed(f.size, 1)}σ"
    if f.kind == "drift":
        return f"dérive de {signed(f.size, 2)}σ par sous-groupe"
    return f"écart-type × {f.size:g}".replace(".", ",")


def fr(x: float, decimals: int = 1) -> str:
    return f"{x:,.{decimals}f}".replace(",", " ").replace(".", ",")


def scenario(rng: np.random.Generator) -> dict[str, dict]:
    """Phase I (sous contrôle) puis phase II (défaut à `start`, compté depuis le début de la phase II)."""
    out = {}
    for key, fault in SCENARIO.items():
        c = CHARACTERISTICS[key]
        phase1 = simulate(c, PHASE1, N, rng)
        phase2 = simulate(c, PHASE2, N, rng, fault)
        out[key] = {"c": c, "fault": fault, "limits": xbar_r_limits(phase1), "phase1": phase1, "phase2": phase2}
    return out


def detection_table(data: dict, params: dict) -> pd.DataFrame:
    """Délai de détection (en sous-groupes après le début du défaut) de chaque méthode."""
    rows = []
    for key, d in data.items():
        lim, fault = d["limits"], d["fault"]
        z = lim.standardize(d["phase2"].mean(axis=1))
        r_flags = np.ptp(d["phase2"], axis=1) > lim.r_ucl
        row = {"caractéristique": d["c"].label, "défaut": fault_label(fault)}
        for m, label in METHODS.items():
            flags = alarms(m, z, **params.get(m, {}))
            false_before = int(flags[:fault.start].sum())
            rank = first_alarm(flags[fault.start:])
            row[label] = "non détecté" if rank > PHASE2 - fault.start else f"{int(rank)}"
            row[f"fausses alarmes {label}"] = false_before
        # La carte R surveille la dispersion : elle n'a de sens que pour le défaut de variance.
        rank_r = first_alarm(r_flags[fault.start:])
        if fault.kind != "variance":
            row["carte R"] = "—"
        else:
            row["carte R"] = "non détecté" if rank_r > PHASE2 - fault.start else f"{int(rank_r)}"
        rows.append(row)
    return pd.DataFrame(rows)


def fig_charts(data: dict, path: Path) -> None:
    """Petits multiples : la carte adaptée à chaque défaut (X̄ pour la moyenne, R pour la dispersion)."""
    fig, axes = plt.subplots(3, 1, figsize=(11, 8.2), sharex=True)
    t = np.arange(1, PHASE1 + PHASE2 + 1)
    for ax, key in zip(axes, ("humidite", "resistance", "granulometrie")):
        d = data[key]
        lim, c, fault = d["limits"], d["c"], d["fault"]
        both = np.vstack([d["phase1"], d["phase2"]])
        if fault.kind == "variance":
            values, center, lcl, ucl, what = np.ptp(both, axis=1), lim.r_center, lim.r_lcl, lim.r_ucl, "Étendue R"
        else:
            values, center, lcl, ucl, what = both.mean(axis=1), lim.center, lim.xbar_lcl, lim.xbar_ucl, "Moyenne X̄"
        out = (values > ucl) | (values < lcl)
        ax.plot(t, values, color=style.SERIES[0], linewidth=1.2, marker="o", markersize=3)
        ax.scatter(t[out], values[out], s=36, color=style.CRITICAL, zorder=3, edgecolors=style.SURFACE,
                   linewidths=0.8, label="Point hors limites")
        for y, name in ((ucl, "LSC"), (center, "LC"), (lcl, "LIC")):
            if name == "LIC" and y <= 0:
                continue
            ax.axhline(y, color=style.INK_2 if name != "LC" else style.BASELINE, linewidth=0.9)
            ax.annotate(name, (t[-1], y), xytext=(6, 0), textcoords="offset points", va="center",
                        fontsize=8, color=style.INK_2, annotation_clip=False)
        onset = PHASE1 + fault.start + 0.5
        ax.axvline(PHASE1 + 0.5, color=style.BASELINE, linewidth=0.8)
        ax.axvline(onset, color=style.INK_2, linewidth=0.8)
        ax.annotate(f"Début du défaut : {fault_label(fault)}", (onset, ax.get_ylim()[1]), xytext=(4, -2),
                    textcoords="offset points", va="top", fontsize=8.5, color=style.INK_2)
        ax.set_title(f"{c.label} ({c.unit}) · {what}", fontsize=10.5)
    axes[0].annotate("Phase I : estimation des limites", (PHASE1 / 2, axes[0].get_ylim()[0]), xytext=(0, 4),
                     textcoords="offset points", ha="center", fontsize=8.5, color=style.MUTED)
    axes[0].legend(loc="upper left")
    axes[-1].set_xlabel("Sous-groupe (4 prélèvements)")
    fig.suptitle("Cartes de contrôle sur un scénario simulé de 130 sous-groupes", x=0.06, ha="left",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_arl(table: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    shifts = table.index.to_numpy()
    for m, label in METHODS.items():
        ax.plot(shifts, table[label].to_numpy(), color=COLORS[m], marker="o", markersize=4, label=label)
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(lambda v, _: fr(v, 2))
    ax.yaxis.set_major_formatter(lambda v, _: fr(v, 0))
    ax.set_xlabel("Décalage de la moyenne du procédé (en écarts-types σ)")
    ax.set_ylabel("ARL : sous-groupes avant alarme")
    ax.set_title("Délai moyen de détection, à fausses alarmes égales pour Shewhart, EWMA et CUSUM")
    ax.annotate("Les règles WE détectent plus vite,\nmais déclenchent 4 fois plus de fausses alarmes",
                (0.0, table.loc[0.0, METHODS["western_electric"]]), xytext=(0.62, 160), textcoords="data",
                fontsize=8.5, color=style.INK_2, va="center",
                arrowprops={"arrowstyle": "-", "color": style.MUTED, "lw": 0.8, "shrinkB": 4})
    ax.legend(loc="upper right")
    fig.savefig(path)
    plt.close(fig)


def fig_detection(params: dict, path: Path) -> pd.DataFrame:
    """Probabilité d'avoir détecté un décalage de 0,5σ après k sous-groupes, sur 4 000 simulations.

    Plus robuste qu'une trajectoire unique, où le hasard peut avantager n'importe quelle méthode.
    """
    horizon = 80
    k = np.arange(1, horizon + 1)
    fig, ax = plt.subplots(figsize=(9, 4.6))
    rows = []
    for m, label in METHODS.items():
        rl = run_lengths(m, DETECTION_SHIFT * np.sqrt(N), n_runs=4000, horizon=horizon, seed=23,
                         **params.get(m, {}))
        ax.plot(k, (rl[:, None] <= k).mean(axis=0), color=COLORS[m], label=label)
        rows.append({"méthode": label, "délai médian": f"{np.median(rl):.0f}",
                     "détecté en 5 sous-groupes": f"{fr(100 * (rl <= 5).mean(), 0)} %",
                     "en 10": f"{fr(100 * (rl <= 10).mean(), 0)} %",
                     "en 20": f"{fr(100 * (rl <= 20).mean(), 0)} %"})
    table = pd.DataFrame(rows)
    ax.axvline(10, color=style.BASELINE, linewidth=0.9)
    ax.annotate(f"Après 10 sous-groupes : EWMA {table.loc[2, 'en 10']}, Shewhart {table.loc[0, 'en 10']}",
                (10, 0.04), xytext=(6, 0), textcoords="offset points", fontsize=8.5, color=style.INK_2)
    ax.set_ylim(0, 1.02)
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:.0%}".replace("%", " %"))
    ax.set_xlabel("Sous-groupes écoulés depuis le début du défaut")
    ax.set_ylabel("Défaut détecté")
    ax.set_title(f"Décalage de {signed(DETECTION_SHIFT, 1)}σ : part des 4 000 simulations où le défaut est détecté")
    ax.legend(loc="lower right")
    fig.savefig(path)
    plt.close(fig)
    return table


def main() -> None:
    style.apply()
    FIGURES.mkdir(parents=True, exist_ok=True)

    # 1. Calibrage : EWMA et CUSUM réglés sur le même ARL0 que Shewhart (≈ 370 sous-groupes).
    L = calibrate("ewma", "L", 2.0, 3.5, lam=LAMBDA)
    h = calibrate("cusum", "h", 3.0, 6.0, k=K)
    params = {"ewma": {"lam": LAMBDA, "L": L}, "cusum": {"k": K, "h": h}}

    # 2. ARL selon le décalage de la moyenne du procédé (n = 4 : décalage réduit = δ·√n).
    shifts = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
    rows, errors = [], []
    for delta in shifts:
        row, err = {}, {}
        for m, label in METHODS.items():
            value, se = arl(m, delta * np.sqrt(N), seed=11, **params.get(m, {}))
            row[label], err[label] = value, se
        rows.append(row)
        errors.append(err)
    arl_table = pd.DataFrame(rows, index=pd.Index(shifts, name="décalage δ (σ)"))
    exact = [shewhart_arl(d * np.sqrt(N)) for d in shifts]

    # 3. Scénario de phase II, délais de détection et capabilité estimée en phase I.
    data = scenario(np.random.default_rng(SEED))
    detection = detection_table(data, params)
    capab = []
    for key, d in data.items():
        c = d["c"]
        cp, cpk = capability(d["phase1"], c.lsl, c.usl, d["limits"].sigma)
        capab.append({"caractéristique": c.label, "tolérance": f"{fr(c.lsl)} – {fr(c.usl)} {c.unit}",
                      "σ estimé (R̄/d2)": fr(d["limits"].sigma, 2), "Cp": fr(cp, 2), "Cpk": fr(cpk, 2)})

    fig_charts(data, FIGURES / "cartes_controle.png")
    fig_arl(arl_table, FIGURES / "arl.png")
    detection_dist = fig_detection(params, FIGURES / "detection.png")

    arl_md = arl_table.copy()
    for col in arl_md:
        arl_md[col] = [f"{fr(v)} ± {fr(1.96 * e[col])}" for v, e in zip(arl_table[col], errors)]
    arl_md.insert(0, "Shewhart exact", [fr(v) for v in exact])
    arl_md.index = [fr(i, 2) for i in arl_md.index]
    det_cols = ["caractéristique", "défaut", *METHODS.values(), "carte R"]
    lines = [
        "# Résultats",
        "",
        f"Calibrage sur ARL0 = {fr(TARGET_ARL0)} (Shewhart 3σ) : EWMA λ = 0,2 → L = {fr(L, 3)} ; "
        f"CUSUM k = 0,5 → h = {fr(h, 3)}. 4 000 séries simulées par point, sous-groupes de taille 4.",
        "",
        "## ARL selon le décalage de la moyenne (intervalle de confiance à 95 %)",
        "",
        arl_md.to_markdown(),
        "",
        f"## Décalage de {signed(DETECTION_SHIFT, 1)}σ : délais de détection sur 4 000 simulations",
        "",
        detection_dist.to_markdown(index=False),
        "",
        "## Scénario simulé : délai de détection (sous-groupes après le début du défaut)",
        "",
        detection[det_cols].to_markdown(index=False),
        "",
        "Fausses alarmes avant le défaut : "
        + "; ".join(f"{d['caractéristique']} : " + ", ".join(
            f"{label} {d[f'fausses alarmes {label}']}" for label in METHODS.values())
            for d in detection.to_dict("records")) + ".",
        "",
        "## Capabilité (phase I)",
        "",
        pd.DataFrame(capab).to_markdown(index=False),
        "",
    ]
    (REPORTS / "resultats.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
