import numpy as np
import pytest

from spcfonderie.arl import TARGET_ARL0, arl, shewhart_arl
from spcfonderie.charts import (alarms, capability, cusum, ewma, ewma_limit, first_alarm,
                                western_electric, xbar_r_limits)
from spcfonderie.process import CHARACTERISTICS, Fault, simulate


def test_xbar_r_limits_recover_process():
    rng = np.random.default_rng(0)
    phase1 = 10 + 2 * rng.standard_normal((2000, 4))
    lim = xbar_r_limits(phase1)
    assert lim.center == pytest.approx(10, abs=0.05)
    assert lim.sigma == pytest.approx(2, rel=0.03)              # R̄ / d2 estime σ sans biais
    assert lim.xbar_ucl - lim.center == pytest.approx(3 * lim.sigma / 2, rel=0.01)   # A2·R̄ = 3σ/√n


@pytest.mark.parametrize("pattern, rule", [
    ([0, 0, 3.5], "1 point au-delà de 3σ"),
    ([0, 2.5, 0.3, 2.2], "2 sur 3 au-delà de 2σ"),
    ([0, 1.2, 1.5, 0.2, 1.1, 1.3], "4 sur 5 au-delà de 1σ"),
    ([0.1] * 8, "8 points du même côté"),
])
def test_western_electric_detects_pattern_at_last_point(pattern, rule):
    flags = western_electric(np.array(pattern))[rule]
    assert flags[-1]
    assert not flags[:-1].any()


def test_western_electric_quiet_on_centered_series():
    z = np.tile([0.5, -0.5], 20)
    assert not np.logical_or.reduce(list(western_electric(z).values())).any()


def test_rules_require_same_side():
    z = np.array([2.5, -2.5, 0.0])
    assert not western_electric(z)["2 sur 3 au-delà de 2σ"].any()


def test_ewma_and_cusum_track_a_shift():
    z = np.r_[np.zeros(20), np.ones(30)]
    assert ewma(z, 0.2)[-1] == pytest.approx(1, abs=0.01)
    hi, lo = cusum(z, 0.5)
    assert hi[-1] == pytest.approx(15)                         # 30 × (1 − 0,5)
    assert lo.max() == 0


def test_first_alarm_rank_and_no_alarm():
    flags = np.array([[False, True, True], [False, False, False]])
    assert first_alarm(flags).tolist() == [2, 4]


def test_shewhart_arl_formula_and_simulation_agree():
    assert shewhart_arl(0) == pytest.approx(370.4, abs=0.1)
    simulated, se = arl("shewhart", 1.0, n_runs=6000, horizon=600)
    assert simulated == pytest.approx(shewhart_arl(1.0), abs=4 * se)


def test_western_electric_rules_cost_false_alarms():
    in_control, _ = arl("western_electric", 0.0, n_runs=3000, horizon=1500)
    assert in_control == pytest.approx(91.75, rel=0.08)        # valeur de référence de la littérature


def test_ewma_limit_formula():
    assert ewma_limit(1.0, 3) == pytest.approx(3)              # λ = 1 : on retrouve Shewhart


def test_capability_indices():
    x = np.full(100, 11.0)
    cp, cpk = capability(x, lsl=4, usl=16, sigma=1)
    assert cp == pytest.approx(2.0)
    assert cpk == pytest.approx(5 / 3)


def test_simulated_faults():
    c = CHARACTERISTICS["humidite"]
    rng = np.random.default_rng(1)
    shifted = simulate(c, 4000, 4, rng, Fault("shift", 2000, 1.0))
    assert shifted[2000:].mean() - shifted[:2000].mean() == pytest.approx(c.sigma, rel=0.05)
    wider = simulate(c, 4000, 4, rng, Fault("variance", 2000, 2.0))
    assert wider[2000:].std() / wider[:2000].std() == pytest.approx(2.0, rel=0.05)
    with pytest.raises(ValueError):
        simulate(c, 10, 4, rng, Fault("inconnu", 0, 1))


def test_alarms_rejects_unknown_method():
    with pytest.raises(ValueError):
        alarms("inconnue", np.zeros(5))


def test_target_arl0():
    assert TARGET_ARL0 == pytest.approx(370.4, abs=0.1)
