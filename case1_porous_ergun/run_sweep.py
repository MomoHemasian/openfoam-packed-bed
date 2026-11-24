#!/usr/bin/env python3
"""Case 1: run porousSimpleFoam for several gas velocities and compare the
pressure drop over the bed with the Ergun equation.

Usage (inside an OpenFOAM terminal):  python3 run_sweep.py
"""
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

# ---- inputs (change these and re-run) ------------------------------------
EPS = 0.40          # bed porosity [-]
DP = 0.012          # pellet diameter [m]
L_BED = 0.50        # bed height [m]  (must match system/topoSetDict)
RHO = 1.204         # air density [kg/m3]
NU = 1.516e-5       # air kinematic viscosity [m2/s]
MU = RHO * NU
VELOCITIES = [0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5]   # superficial, m/s
# --------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"


def ergun_dp(u):
    """Ergun pressure drop [Pa] over the bed for superficial velocity u."""
    viscous = 150 * MU * (1 - EPS) ** 2 / (EPS ** 3 * DP ** 2) * u
    inertial = 1.75 * RHO * (1 - EPS) / (EPS ** 3 * DP) * u ** 2
    return (viscous + inertial) * L_BED


def darcy_forchheimer_coeffs():
    d = 150 * (1 - EPS) ** 2 / (EPS ** 3 * DP ** 2)   # 1/m2
    f = 2 * 1.75 * (1 - EPS) / (EPS ** 3 * DP)        # 1/m
    return d, f


def sh(cmd, cwd):
    with open(cwd / f"log.{cmd.split()[0]}", "w") as log:
        subprocess.run(cmd.split(), cwd=cwd, stdout=log, stderr=subprocess.STDOUT, check=True)


def read_internal_field(path):
    """Read the cell values of a scalar OpenFOAM field written in ascii."""
    text = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)\s*;", text, re.S)
    return np.array(m.group(2).split(), dtype=float)


def read_column_p(case):
    """The mesh is 400 equal cells along z (0 to 0.8 m), so cell i sits at z = (i + 0.5) * 0.002."""
    times = [d for d in case.iterdir() if d.is_dir() and d.name.replace(".", "").isdigit() and d.name != "0"]
    last = max(times, key=lambda d: float(d.name))
    p = read_internal_field(last / "p")
    z = (np.arange(p.size) + 0.5) * 0.8 / p.size
    return z, p


def run_one(u):
    case = RUNS / f"U_{u:.2f}"
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    for sub in ("0", "constant", "system"):
        shutil.copytree(HERE / sub, case / sub)
    d, f = darcy_forchheimer_coeffs()
    tmpl = case / "constant" / "porosityProperties.template"
    text = tmpl.read_text().replace("D_VALUE", f"{d:.6e}").replace("F_VALUE", f"{f:.6e}")
    (case / "constant" / "porosityProperties").write_text(text)
    tmpl.unlink()
    for name in ("U",):
        p = case / "0" / name
        p.write_text(p.read_text().replace("U_IN", f"{u}"))
    sh("blockMesh", case)
    sh("topoSet", case)
    sh("porousSimpleFoam", case)
    # read the cell values of p at the last time and take p(z=0.1) - p(z=0.6)
    z, p_kin = read_column_p(case)
    # Inside the bed p falls linearly. Fit a straight line to the cells well inside
    # the bed and multiply the slope by the bed height (avoids the kink at the bed ends).
    inside = (z > 0.12) & (z < 0.58)
    slope = np.polyfit(z[inside], p_kin[inside], 1)[0]
    dp_bed = -slope * L_BED * RHO
    return dp_bed, z, p_kin * RHO


def main():
    d, f = darcy_forchheimer_coeffs()
    print(f"eps = {EPS}, dp = {DP} m  ->  d = {d:.4e} 1/m2, f = {f:.4e} 1/m")
    rows = []
    profiles = {}
    for u in VELOCITIES:
        dp_cfd, z, p = run_one(u)
        dp_erg = ergun_dp(u)
        re_p = RHO * u * DP / (MU * (1 - EPS))
        err = 100 * (dp_cfd - dp_erg) / dp_erg
        rows.append((u, re_p, dp_cfd, dp_erg, err))
        profiles[u] = (z, p)
        print(f"U = {u:5.2f} m/s  Re_p = {re_p:7.1f}  dp OpenFOAM = {dp_cfd:9.2f} Pa  "
              f"dp Ergun = {dp_erg:9.2f} Pa  diff = {err:+.3f} %")
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    np.savetxt(out / "pressure_drop.csv", np.array(rows), delimiter=",",
               header="U_superficial_m_s,Re_p_modified,dp_OpenFOAM_Pa,dp_Ergun_Pa,difference_percent",
               comments="", fmt="%.6g")
    np.savez(out / "profiles.npz", **{f"U_{u:.2f}": np.vstack(v) for u, v in profiles.items()})
    print("Results written to", out)


if __name__ == "__main__":
    main()
