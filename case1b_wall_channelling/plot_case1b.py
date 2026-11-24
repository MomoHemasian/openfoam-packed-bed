#!/usr/bin/env python3
"""Plot for Case 1b: gas velocity across a 2D grate bed with a looser strip at the side wall."""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "plots"
NX, NZ, W, H = 200, 110, 0.5, 0.55          # mesh from system/blockMeshDict
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"


def read_u(case):
    times = [d for d in case.iterdir() if d.is_dir() and d.name.isdigit() and d.name != "0"]
    text = (max(times, key=lambda d: int(d.name)) / "U").read_text()
    body = re.search(r"internalField\s+nonuniform\s+List<vector>\s*\d+\s*\((.*?)\n\)\s*;", text, re.S).group(1)
    vec = np.array([v.split() for v in re.findall(r"\(([^()]+)\)", body)], dtype=float)
    return vec.reshape(NZ, NX, 3)          # blockMesh numbers cells x first, then z


def main():
    U = read_u(HERE)
    x = (np.arange(NX) + 0.5) * W / NX
    z = (np.arange(NZ) + 0.5) * H / NZ
    uz = U[:, :, 2]
    row = uz[np.argmin(abs(z - 0.40))]
    core = row[x < 0.40].mean()
    strip = row[x > 0.48].mean()
    share = row[x > 0.48].sum() / row.sum()
    print(f"core u = {core:.3f} m/s, wall strip u = {strip:.3f} m/s, ratio = {strip / core:.2f}, "
          f"strip carries {100 * share:.1f} % of the gas on {100 * 0.02 / W:.0f} % of the width")

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4), dpi=150, gridspec_kw={"width_ratios": [1.1, 1]})
    im = a1.pcolormesh(x, z, uz, cmap="Blues", shading="auto")
    for zz in (0.05, 0.15):
        a1.axhline(zz, color=MUTED, lw=0.8, ls="--")
    a1.axvline(0.48, color=MUTED, lw=0.8, ls="--")
    a1.text(0.02, 0.10, "hearth layer", color=INK, va="center")
    a1.text(0.02, 0.35, "green pellets", color=INK, va="center")
    a1.text(0.02, 0.025, "wind box", color=INK, va="center")
    a1.set_xlabel("Distance from grate centre x (m)")
    a1.set_ylabel("Height z (m)")
    a1.set_title("Vertical gas velocity (m/s)", fontsize=9.5, loc="left", color=INK)
    fig.colorbar(im, ax=a1, fraction=0.05)
    a2.plot(x, row, color=BLUE, lw=2)
    a2.axvspan(0.48, 0.5, color="#fbd9cb", lw=0)
    a2.text(0.475, core * 1.05, "looser strip\nat the wall", ha="right", color=MUTED, fontsize=9)
    a2.set_xlabel("Distance from grate centre x (m)")
    a2.set_ylabel("Vertical gas velocity at z = 0.40 m (m/s)")
    a2.set_ylim(0, row.max() * 1.15)
    a2.grid(True, color="#e4e3df", lw=0.6)
    a2.set_title(f"Gas velocity is {100 * (strip / core - 1):.0f} % higher in the wall strip",
                 fontsize=9.5, loc="left", color=INK)
    fig.suptitle("Case 1b: 2D grate bed, 2000 Pa across the bed, air 20 C", x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "case1b_wall_channelling.png")


if __name__ == "__main__":
    main()
