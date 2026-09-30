"""SynchronyLab — reference Python implementation.

A two-region (hippocampus, prefrontal cortex) spiking network of Izhikevich
neurons with pooled synapses, plus the spectral analyses used in the browser
app: Welch power spectra, magnitude-squared coherence and the Tort modulation
index for theta–gamma phase–amplitude coupling.

The model mirrors ``index.html`` in structure and parameters. Random number
streams differ between the two, so individual runs are not bit-identical, but
the statistics are the same.

Author: Mersedeh Shamsedin — MIT License
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import coherence, welch

NE, NI = 160, 40          # pyramidal (RS) and fast-spiking cells per region
DT = 0.5                  # ms
FS = 1000.0               # LFP sampling rate (Hz)
TAU_E, TAU_I = 2.0, 5.0   # synaptic decay (ms)
THETA_HZ = 8.0            # septal drive frequency
THETA_BAND = (6.0, 10.0)
GAMMA_BAND = (30.0, 80.0)


@dataclass
class Params:
    integrity: float = 1.0   # 1 = intact synapses; lower = synapse loss
    coupling: float = 1.0    # long-range HPC <-> PFC gain
    theta: float = 1.0       # septal theta drive gain
    drive: float = 1.0       # background excitation gain
    noise: float = 1.0       # membrane noise gain


@dataclass
class Region:
    rng: np.random.Generator
    a: np.ndarray = field(init=False)
    b: np.ndarray = field(init=False)
    c: np.ndarray = field(init=False)
    d: np.ndarray = field(init=False)
    v: np.ndarray = field(init=False)
    u: np.ndarray = field(init=False)
    k: np.ndarray = field(init=False)
    bias: np.ndarray = field(init=False)
    sE: float = 0.0
    sI: float = 0.0

    def __post_init__(self) -> None:
        r = self.rng.random(NE + NI)
        re, ri = r[:NE], r[NE:]
        self.a = np.r_[np.full(NE, 0.02), 0.02 + 0.08 * ri]
        self.b = np.r_[np.full(NE, 0.2), 0.25 - 0.05 * ri]
        self.c = np.r_[-65 + 15 * re**2, np.full(NI, -65.0)]
        self.d = np.r_[8 - 6 * re**2, np.full(NI, 2.0)]
        self.v = np.full(NE + NI, -65.0)
        self.u = self.b * self.v
        self.k = 0.7 + 0.6 * self.rng.random(NE + NI)
        self.bias = self.rng.standard_normal(NE + NI)
        self.exc = np.arange(NE + NI) < NE


class Network:
    """Hippocampus → prefrontal spiking network."""

    def __init__(self, params: Params | None = None, seed: int = 0) -> None:
        self.p = params or Params()
        self.rng = np.random.default_rng(seed)
        self.H = Region(self.rng)
        self.F = Region(self.rng)
        self.t = 0.0

    def _step_region(self, R: Region, ext: float, ext_i: float,
                     theta_e: float, theta_i: float):
        p = self.p
        q = p.integrity
        wEE, wEI, wIE, wII = 6 * q, 14 * (0.35 + 0.65 * q), 9.0, 4.0
        syn_e = wEE * R.sE - wIE * R.sI
        syn_i = wEI * R.sE - wII * R.sI
        base = np.where(R.exc, 4.6, 3.0) * p.drive
        I = (base + 0.8 * R.bias
             + p.noise * 2.2 * self.rng.standard_normal(NE + NI)
             + np.where(R.exc, R.k * syn_e + ext + theta_e,
                        R.k * syn_i + ext_i + theta_i))
        v, u = R.v, R.u
        for _ in range(2):  # two half-steps for v (Izhikevich 2003)
            v = v + 0.5 * DT * (0.04 * v * v + 5 * v + 140 - u + I)
        u = u + DT * R.a * (R.b * v - u)
        fired = v >= 30
        R.v = np.where(fired, R.c, v)
        R.u = np.where(fired, u + R.d, u)
        lfp = syn_e + ext
        return fired, lfp

    def step(self):
        p, H, F = self.p, self.H, self.F
        th = np.sin(2 * np.pi * THETA_HZ * self.t / 1000)
        c_hf = 8 * p.coupling * p.integrity
        c_fh = 1.5 * p.coupling * p.integrity
        ext_h, ext_f = c_fh * F.sE, c_hf * H.sE
        fh, lh = self._step_region(H, ext_h, 0.3 * ext_h,
                                   1.6 * p.theta * th, 3.0 * p.theta * th)
        ff, lf = self._step_region(F, ext_f, 0.6 * ext_f, 0.0, 0.0)
        for R, f in ((H, fh), (F, ff)):
            R.sE += -R.sE * DT / TAU_E + f[:NE].sum() * 20 / NE
            R.sI += -R.sI * DT / TAU_I + f[NE:].sum() * 20 / NI
        self.t += DT
        return fh, ff, lh, lf

    def run(self, duration_ms: float, warmup_ms: float = 500.0):
        """Simulate and return LFP proxies (1 kHz) and spike times."""
        for _ in range(int(warmup_ms / DT)):
            self.step()
        n = int(duration_ms / DT)
        lfp_h = np.zeros(n // 2)
        lfp_f = np.zeros(n // 2)
        spikes_h, spikes_f = [], []
        for s in range(n):
            fh, ff, lh, lf = self.step()
            lfp_h[s // 2] += lh / 2
            lfp_f[s // 2] += lf / 2
            t = s * DT
            spikes_h.extend((t, i) for i in np.flatnonzero(fh))
            spikes_f.extend((t, i) for i in np.flatnonzero(ff))
        return {"lfp_h": lfp_h, "lfp_f": lfp_f,
                "spikes_h": np.array(spikes_h), "spikes_f": np.array(spikes_f),
                "duration_s": duration_ms / 1000}


# ---------------------------------------------------------------- analysis
def _band_mean(f, x, band):
    m = (f >= band[0]) & (f <= band[1])
    return float(x[m].mean())


def _analytic_band(x, band):
    """Ideal FFT band-pass followed by Hilbert transform."""
    x = x - x.mean()
    X = np.fft.fft(x)
    f = np.fft.fftfreq(len(x), 1 / FS)
    keep = (f > 0) & (f >= band[0]) & (f <= band[1])
    return np.fft.ifft(np.where(keep, 2 * X, 0))


def modulation_index(phase_sig, amp_sig, n_bins=18, edge=300):
    """Tort et al. (2010) modulation index."""
    ph = np.angle(_analytic_band(phase_sig, THETA_BAND))[edge:-edge]
    amp = np.abs(_analytic_band(amp_sig, GAMMA_BAND))[edge:-edge]
    bins = np.minimum(((ph + np.pi) / (2 * np.pi) * n_bins).astype(int), n_bins - 1)
    mean_amp = np.array([amp[bins == b].mean() for b in range(n_bins)])
    p = mean_amp / mean_amp.sum()
    H = -(p * np.log(p)).sum()
    return float((np.log(n_bins) - H) / np.log(n_bins)), p


def analyse(lfp_h, lfp_f):
    f, pxx = welch(lfp_h, FS, window="hann", nperseg=512, noverlap=256)
    _, pyy = welch(lfp_f, FS, window="hann", nperseg=512, noverlap=256)
    _, coh = coherence(lfp_h, lfp_f, FS, window="hann", nperseg=512, noverlap=256)
    g = (f >= 25) & (f <= 90)
    mi, _ = modulation_index(lfp_h, lfp_h)
    mi_x, _ = modulation_index(lfp_h, lfp_f)
    return {
        "f": f, "psd_h": pxx, "psd_f": pyy, "coherence": coh,
        "theta_coherence": _band_mean(f, coh, THETA_BAND),
        "gamma_coherence": _band_mean(f, coh, GAMMA_BAND),
        "gamma_peak_hz": float(f[g][np.argmax(pxx[g])]),
        "pac_mi_hpc": mi, "pac_mi_hpc_to_pfc": mi_x,
    }
