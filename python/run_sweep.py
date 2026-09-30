"""Reproduce the synaptic-integrity sweep and save a publication-style figure.

Usage:
    python python/run_sweep.py --seeds 5 --duration 8 --out docs/sweep.png
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from synchronylab import Network, Params, analyse  # noqa: E402

LEVELS = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--duration", type=float, default=8.0, help="seconds per run")
    ap.add_argument("--coupling", type=float, default=1.0)
    ap.add_argument("--theta", type=float, default=1.0)
    ap.add_argument("--out", default="docs/sweep.png")
    args = ap.parse_args()

    rows = []
    for q in LEVELS:
        for s in range(args.seeds):
            net = Network(Params(integrity=q, coupling=args.coupling, theta=args.theta), seed=1000 + s)
            r = analyse(**{k: v for k, v in net.run(args.duration * 1000).items() if k.startswith("lfp")})
            rows.append({"integrity": q, "seed": s, **{k: r[k] for k in
                         ("theta_coherence", "gamma_coherence", "pac_mi_hpc", "pac_mi_hpc_to_pfc", "gamma_peak_hz")}})
            print(f"integrity={q:.1f} seed={s} θcoh={r['theta_coherence']:.2f} "
                  f"γcoh={r['gamma_coherence']:.2f} MI={r['pac_mi_hpc']*1e3:.1f}e-3", flush=True)

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out.with_suffix(".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def agg(key):
        m = np.array([np.mean([r[key] for r in rows if r["integrity"] == q]) for q in LEVELS])
        se = np.array([np.std([r[key] for r in rows if r["integrity"] == q], ddof=1)
                       / np.sqrt(args.seeds) if args.seeds > 1 else 0 for q in LEVELS])
        return m, se

    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    for key, lab, col in (("theta_coherence", "θ coherence (6–10 Hz)", "#c7770f"),
                          ("gamma_coherence", "γ coherence (30–80 Hz)", "#1f79b5")):
        m, se = agg(key)
        axs[0].errorbar(LEVELS, m, yerr=se, marker="o", capsize=3, color=col, label=lab)
    axs[0].set(xlabel="Synaptic integrity", ylabel="HPC–PFC coherence", ylim=(0, 1))
    axs[0].invert_xaxis()
    axs[0].legend(frameon=False)
    m, se = agg("pac_mi_hpc")
    axs[1].errorbar(LEVELS, m * 1e3, yerr=se * 1e3, marker="o", capsize=3, color="#7a52c2")
    axs[1].set(xlabel="Synaptic integrity", ylabel="θ–γ modulation index (×10⁻³)")
    axs[1].invert_xaxis()
    for ax in axs:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"SynchronyLab integrity sweep · mean ± SEM, n = {args.seeds} networks per level")
    fig.savefig(out, dpi=200)
    print("saved", out, "and", out.with_suffix(".csv"))


if __name__ == "__main__":
    main()
