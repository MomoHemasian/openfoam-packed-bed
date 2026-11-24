#!/usr/bin/env python3
"""Plots for Case 1 (porous zone versus Ergun)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from run_sweep import EPS, DP, L_BED, ergun_dp

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "plots"
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e3df",
                     "grid.linewidth": 0.6})


def main():
    OUT.mkdir(exist_ok=True)
    data = np.loadtxt(HERE / "results" / "pressure_drop.csv", delimiter=",", skiprows=1)
    u, dp_cfd = data[:, 0], data[:, 2]

    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    uu = np.linspace(0, u.max() * 1.05, 200)
    ax.plot(uu, ergun_dp(uu) / L_BED / 1000, color=BLUE, lw=2, label="Ergun equation")
    ax.plot(u, dp_cfd / L_BED / 1000, "o", ms=8, color=ORANGE, mec="white", mew=1.5,
            label="OpenFOAM, porous zone")
    ax.set_xlabel("Superficial gas velocity u (m/s)")
    ax.set_ylabel("Pressure gradient in the bed (kPa/m)")
    ax.set_title(f"Case 1: pellet bed as a porous zone (eps = {EPS}, d_p = {DP*1000:.0f} mm, air 20 C)",
                 fontsize=9.5, color=INK, loc="left")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "case1_pressure_drop_vs_ergun.png")

    prof = np.load(HERE / "results" / "profiles.npz")
    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    for key, col in (("U_0.50", BLUE), ("U_1.00", ORANGE)):
        z, p = prof[key]
        ax.plot(z, p, color=col, lw=2, label=f"u = {float(key[2:]):.1f} m/s")
    ax.axvspan(0.1, 0.6, color="#cde2fb", alpha=0.5, lw=0)
    ax.text(0.35, ax.get_ylim()[1] * 0.92, "pellet bed", ha="center", color=MUTED)
    ax.set_xlabel("Height in the column z (m)")
    ax.set_ylabel("Pressure above outlet (Pa)")
    ax.set_title("Case 1: pressure along the column", fontsize=9.5, color=INK, loc="left")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "case1_pressure_profile.png")


if __name__ == "__main__":
    main()
