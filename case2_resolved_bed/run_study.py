#!/usr/bin/env python3
"""Case 2: run the resolved-particle study and post-process it.

For every combination of bed length (in lattice periods), gas velocity and
mesh level it: writes the case (make_case.py), meshes it (blockMesh +
snappyHexMesh), solves the flow (simpleFoam), then solves the temperature in
the frozen flow (scalarTransportFoam), and reads back
  * the pressure at the inlet (the outlet is at p = 0), and
  * the flow-weighted (mixing-cup) temperature at the outlet.

Entrance and exit effects are removed by running two bed lengths, N1 and N2
periods, and using only the DIFFERENCE between them:
  dp/L         = (dp(N2) - dp(N1)) / ((N2 - N1) s)
  ln(T1/T2)    = h a (N2 - N1) s / (rho cp u)      ->  h, then Nu = h d / k

Usage (inside an OpenFOAM terminal):  python3 run_study.py
Runs 4 cases at a time; the whole study takes roughly 30-60 min on 4 cores.
"""
import json
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import make_case as mc

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
RE_LIST = [1, 5, 20, 50, 100]          # particle Reynolds number rho u d / mu (superficial u)
PERIODS = (2, 4)
LEVELS = (2, 3)                        # mesh refinement at the pellet surface (3 = finer)
LEVEL3_RE = (5, 50)                    # the fine mesh is run only at these Re (mesh check)


def sh(cmd, cwd):
    name = cmd.split()[0]
    with open(cwd / f"log.{name}", "w") as log:
        subprocess.run(cmd.split(), cwd=cwd, stdout=log, stderr=subprocess.STDOUT, check=True)


def last_time_dir(case, field=None):
    times = [d for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.eE+-]+", d.name) and d.name != "0"
             and (field is None or (d / field).exists())]
    return max(times, key=lambda d: float(d.name))


def patch_values(field_file, patch):
    """Return the face values of one patch from an ascii OpenFOAM field file."""
    text = field_file.read_text()
    block = re.search(patch + r"\s*\{(.*?)\}", text.split("boundaryField", 1)[1], re.S).group(1)
    m = re.search(r"value\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)\s*;", block, re.S)
    if m:
        return np.array(m.group(1).split(), dtype=float)
    m = re.search(r"value\s+uniform\s+([-0-9.eE+]+)\s*;", block)
    return np.array([float(m.group(1))])


def run_one(args):
    periods, re_p, level = args
    u = re_p * mc.NU / mc.D
    case = RUNS / f"L{level}_N{periods}_Re{re_p}"
    done = case / "result.json"
    if done.exists():
        return json.loads(done.read_text())
    finished = (case / "log.scalarTransportFoam").exists() and \
        "End" in (case / "log.scalarTransportFoam").read_text()
    if not finished:
        mc.make(case, periods, u, level, base_cells=10)
        sh("blockMesh", case)
        sh("snappyHexMesh -overwrite", case)
        sh("checkMesh", case)
        sh("simpleFoam", case)
        t = last_time_dir(case)
        (t / "T").write_text((case / "0" / "T").read_text())   # start T from the converged flow
        sh("scalarTransportFoam", case)
    t_flow = last_time_dir(case, "p")                       # simpleFoam result
    t_temp = last_time_dir(case, "T")                       # scalarTransportFoam result
    p_in = patch_values(t_flow / "p", "inlet").mean()       # inlet faces are all the same size
    phi = patch_values(t_flow / "phi", "outlet")
    temp = patch_values(t_temp / "T", "outlet")
    t_mix = float((phi * temp).sum() / phi.sum())
    n_cells = int(re.search(r"cells:\s+(\d+)", (case / "log.checkMesh").read_text()).group(1))
    iters = int(re.findall(r"^Time = (\d+)", (case / "log.simpleFoam").read_text(), re.M)[-1])
    converged = "converged" in (case / "log.simpleFoam").read_text()
    res = dict(periods=periods, Re=re_p, level=level, u=u, p_in_kin=float(p_in), T_out=t_mix,
               cells=n_cells, iterations=iters, converged=converged)
    done.write_text(json.dumps(res, indent=1))
    return res


def main():
    jobs = [(n, r, 2) for r in RE_LIST for n in PERIODS]
    jobs += [(n, r, 3) for r in LEVEL3_RE for n in PERIODS]
    RUNS.mkdir(exist_ok=True)
    with ProcessPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(run_one, jobs))
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / "raw_runs.json").write_text(json.dumps(results, indent=1))
    for r in results:
        print(r)


if __name__ == "__main__":
    main()
