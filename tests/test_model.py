import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "python"))
from synchronylab import FS, NE, Network, Params, analyse, modulation_index  # noqa: E402


def test_firing_rates_are_physiological():
    out = Network(Params(), seed=1).run(2000)
    rate_h = (out["spikes_h"][:, 1] < NE).sum() / NE / out["duration_s"]
    assert 2 < rate_h < 30


def test_gamma_emerges_in_intact_network():
    out = Network(Params(), seed=2).run(4096)
    r = analyse(out["lfp_h"], out["lfp_f"])
    assert 25 <= r["gamma_peak_hz"] <= 60


def test_modulation_index_detects_coupling():
    t = np.arange(8192) / FS
    theta = np.sin(2 * np.pi * 8 * t)
    coupled = theta + (1 + theta) * 0.3 * np.sin(2 * np.pi * 50 * t)
    flat = theta + 0.3 * np.sin(2 * np.pi * 50 * t)
    assert modulation_index(coupled, coupled)[0] > 10 * modulation_index(flat, flat)[0]


def test_synapse_loss_reduces_coherence():
    def coh(q):
        vals = []
        for s in range(3):
            out = Network(Params(integrity=q), seed=10 + s).run(6000)
            vals.append(analyse(out["lfp_h"], out["lfp_f"])["gamma_coherence"])
        return np.mean(vals)
    assert coh(1.0) > coh(0.2) + 0.1
