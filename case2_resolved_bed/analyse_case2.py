#!/usr/bin/env python3
"""Case 2: turn the raw runs into pressure-drop and heat-transfer results,
compare them with Ergun, Wakao-Kaguei and Gunn, and make the plots."""
import glob
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import make_case as mc

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "plots"
RHO, K_GAS, CP, PR = 1.204, 0.0257, 1006.0, 0.71      # air at 20 C
MU = RHO * mc.NU
EPS, D, S = mc.porosity(), mc.D, mc.S
A_SPEC = 6 * (1 - EPS) / D                           # pellet surface per bed volume [1/m]
BLUE, ORANGE, AQUA, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"


def ergun_gradient(u):
    return 150 * MU * (1 - EPS) ** 2 / (EPS ** 3 * D ** 2) * u + 1.75 * RHO * (1 - EPS) / (EPS ** 3 * D) * u ** 2


def wakao_kaguei(re):
    return 2 + 1.1 * PR ** (1 / 3) * re ** 0.6


def gunn(re):
    return ((7 - 10 * EPS + 5 * EPS ** 2) * (1 + 0.7 * re ** 0.2 * PR ** (1 / 3))
            + (1.33 - 2.4 * EPS + 1.2 * EPS ** 2) * re ** 0.7 * PR ** (1 / 3))


def load():
    runs = {}
    for f in glob.glob(str(HERE / "runs" / "*" / "result.json")):
        r = json.loads(Path(f).read_text())
        runs[(r["level"], r["periods"], r["Re"])] = r
    return runs


def evaluate(runs, level, re):
    """Use the difference between the 4-period and the 2-period bed (removes entry and exit effects)."""
    a, b = runs[(level, 2, re)], runs[(level, 4, re)]
    u = a["u"]
    length = 2 * S
    grad = (b["p_in_kin"] - a["p_in_kin"]) * RHO / length
    h = RHO * CP * u / (A_SPEC * length) * math.log(a["T_out"] / b["T_out"])
    return dict(Re=re, level=level, u=u, cells_4p=b["cells"], grad_cfd=grad, grad_ergun=ergun_gradient(u),
                ratio=grad / ergun_gradient(u), Nu=h * D / K_GAS, Nu_WK=wakao_kaguei(re), Nu_Gunn=gunn(re),
                T_out_2p=a["T_out"], T_out_4p=b["T_out"], converged=a["converged"] and b["converged"])


def main():
    runs = load()
    rows = []
    for (lvl, n, re) in sorted(runs):
        if n == 2 and (lvl, 4, re) in runs:
            rows.append(evaluate(runs, lvl, re))
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    keys = list(rows[0].keys())
    with open(out / "case2_results.csv", "w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    print(f"porosity = {EPS:.4f}, s = {S * 1000:.3f} mm, a = {A_SPEC:.1f} 1/m")
    for r in rows:
        print(f"level {r['level']}  Re {r['Re']:>4}  cells(4p) {r['cells_4p']:>7}  dp/L CFD {r['grad_cfd']:9.4g} Pa/m  "
              f"Ergun {r['grad_ergun']:9.4g}  ratio {r['ratio']:.3f}  Nu {r['Nu']:6.2f}  WK {r['Nu_WK']:5.2f}  "
              f"Gunn {r['Nu_Gunn']:5.2f}  T_out {r['T_out_2p']:.2e}/{r['T_out_4p']:.2e}  conv {r['converged']}")

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                         "grid.color": "#e4e3df", "grid.linewidth": 0.6})
    l2 = [r for r in rows if r["level"] == 2]
    l3 = [r for r in rows if r["level"] == 3]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4), dpi=150)
    re_line = np.logspace(0, 2, 100)
    u_line = re_line * MU / (RHO * D)
    # dimensionless friction factor: f = (dp/L) d eps^3 / (rho u^2 (1 - eps))
    ff = lambda g, u: g * D * EPS ** 3 / (RHO * u ** 2 * (1 - EPS))
    a1.loglog(re_line, ff(ergun_gradient(u_line), u_line), color=BLUE, lw=2, label="Ergun equation (random beds)")
    a1.loglog([r["Re"] for r in l2], [ff(r["grad_cfd"], r["u"]) for r in l2], "o", ms=8, color=ORANGE,
              mec="white", mew=1.5, label="OpenFOAM, resolved BCC bed")
    if l3:
        a1.loglog([r["Re"] for r in l3], [ff(r["grad_cfd"], r["u"]) for r in l3], "s", ms=7, mfc="none",
                  mec=INK, mew=1.2, label="same, finer mesh")
    a1.set_xlabel("Particle Reynolds number  Re = rho u d / mu")
    a1.set_ylabel("Friction factor  (dp/L) d eps^3 / (rho u^2 (1-eps))")
    a1.set_title("Pressure drop", fontsize=9.5, loc="left", color=INK)
    a1.legend(frameon=False, fontsize=8.5)

    heat = [r for r in l2 if r["Re"] >= 20]
    a2.plot(re_line, wakao_kaguei(re_line), color=BLUE, lw=2, label="Wakao-Kaguei (1982)")
    a2.plot(re_line, gunn(re_line), color=AQUA, lw=2, label="Gunn (1978)")
    a2.plot([r["Re"] for r in heat], [r["Nu"] for r in heat], "o", ms=8, color=ORANGE, mec="white", mew=1.5,
            label="OpenFOAM, resolved BCC bed")
    heat3 = [r for r in l3 if r["Re"] >= 20]
    if heat3:
        a2.plot([r["Re"] for r in heat3], [r["Nu"] for r in heat3], "s", ms=7, mfc="none", mec=INK, mew=1.2,
                label="same, finer mesh")
    a2.set_xscale("log")
    a2.set_xlim(1, 120)
    a2.set_xlabel("Particle Reynolds number  Re = rho u d / mu")
    a2.set_ylabel("Nusselt number  Nu = h d / k")
    a2.set_title("Gas-to-pellet heat transfer (Re >= 20 only, see README)", fontsize=9.5, loc="left", color=INK)
    a2.legend(frameon=False, fontsize=8.5)
    fig.suptitle(f"Case 2: resolved spheres, BCC packing, eps = {EPS:.3f}, d = {D * 1000:.0f} mm, air 20 C",
                 x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "case2_pressure_drop_and_heat_transfer.png")


if __name__ == "__main__":
    main()
