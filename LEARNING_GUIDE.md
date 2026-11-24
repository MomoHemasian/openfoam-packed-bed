# Learning guide: how this project works, step by step

This guide is for me (and for anyone who reads the repository). The goal is
that I can explain and defend every choice in an interview, not only show
plots. Read it next to the case files.

## 0. Installing OpenFOAM on Windows (WSL)

1. In PowerShell (as administrator): `wsl --install -d Ubuntu-24.04`, then restart.
2. In the Ubuntu terminal:
   ```
   sudo apt update
   sudo apt install openfoam python3-numpy python3-matplotlib
   source /usr/share/openfoam/etc/bashrc      # put this line in ~/.bashrc
   ```
   This installs OpenFOAM v1912 (ESI version), the same version used for every
   result in this repository.
3. Copy the repository into the Linux home folder (not into OneDrive, it is
   slow there), for example `cp -r /mnt/c/Users/<me>/Desktop/openfoam-packed-bed ~/`.
4. Run `./Allrun`, or each case on its own (see the README).
5. To look at results in 3D: run `foamToVTK` inside a case and open the `VTK`
   folder in ParaView for Windows (free download).

Known problem of this Ubuntu package: OpenFOAM "function objects" (the
`functions { }` block and `postProcess -func`) stop with an error about
"sha1". That is why all post-processing here is done by reading the
OpenFOAM result files directly in Python. On newer OpenFOAM versions the
function objects work.

## 1. What an OpenFOAM case is

Every case is a folder with three sub-folders:

| Folder | What is inside | My analogy |
|---|---|---|
| `0/` | initial and boundary conditions of each field (`U`, `p`, `T`) | the boundary conditions of my PDEs |
| `constant/` | mesh and physical properties (viscosity, porous coefficients) | the parameters of the model |
| `system/` | mesh recipe, numerical schemes, solvers, run control | the numerical method |

The solver (here `porousSimpleFoam`, `simpleFoam`, `scalarTransportFoam`) is
a program that reads these folders and writes result folders named by
iteration number.

Main difference from my own models: in my induration codes I wrote a
finite-difference scheme with RK4 in time. OpenFOAM uses the finite-volume
method (integrate each equation over each cell, fluxes through faces), and
for steady incompressible flow the SIMPLE algorithm (guess pressure, solve
momentum, correct pressure so that mass is conserved, repeat). My GrateDesigner
solver is already finite-volume, so the idea is the same; what is new is the
pressure-velocity coupling and the unstructured 3D mesh.

## 2. Case 1: the bed as a porous zone

**Idea.** We do not draw the pellets. The bed is a region where the momentum
equation gets an extra sink term (a resistance):

    S = -( mu * d + 0.5 * rho * |U| * f ) * U        (Darcy-Forchheimer)

If I compare this with Ergun,

    dp/L = 150 mu (1-eps)^2 / (eps^3 dp^2) * u + 1.75 rho (1-eps) / (eps^3 dp) * u^2

I get `d = 150 (1-eps)^2 / (eps^3 dp^2)` and `f = 3.5 (1-eps) / (eps^3 dp)`.
`run_sweep.py` computes them from eps and dp.

**Things I must be able to say:**
- `U` in OpenFOAM's porous solver is the *superficial* velocity (flow rate
  divided by the full area), the same velocity that Ergun uses.
- The flow is set to laminar on purpose. In a porous-zone model we only
  resolve the mean flow. The pore-scale inertia and turbulence are already
  inside the Forchheimer (1.75) term.
- Incompressible OpenFOAM solvers use kinematic pressure p/rho, so I multiply
  by rho to get Pascal.
- The result matches Ergun exactly. This is a **verification**: it proves I
  set up the model and the coefficients correctly. It is **not** a
  validation, because I put Ergun in and got Ergun out. The real value of a
  porous model comes in 2D or 3D, where the bed is not uniform (Case 1b).
- In ANSYS Fluent the same model is the "porous zone" with viscous
  resistance `1/alpha` and inertial resistance `C2`. The numbers are exactly
  `1/alpha = d` and `C2 = f`. So this case transfers directly to Fluent.

**Numerics I chose and why.** Under-relaxation 0.3 for p and 0.7 for U
(standard SIMPLE values). My first try with SIMPLEC and no pressure
relaxation diverged, because the porous sink is very strong compared with
the convection term. Relaxation fixed it. Convergence was checked with the
residuals (p below 1e-8).

## 3. Case 1b: why a 2D porous model is useful (wall channelling)

A 2D slice of half a grate (0.5 m wide), with a wind box, a 10 cm hearth
layer and 40 cm of green pellets. Next to the side wall, spheres cannot pack
as tightly, so a 2 cm strip has a higher porosity. The fan pressure is fixed
(2000 Pa across the bed).

Result: the gas is about 40 % faster in the wall strip. This is the known
"wall effect" or edge channelling. On a real induration machine it means the
pellets at the side walls see more gas and heat up differently.

**Honest limits:** the strip porosity (0.48) and width (2 cm) are assumed
values for illustration, not measured. Real near-wall porosity oscillates
over about one to two pellet diameters. The gas is cold air. A real case
needs hot gas (temperature-dependent density and viscosity) and the
temperature field.

## 4. Case 2: the pellets drawn one by one (resolved simulation)

**Geometry.** Spheres of 12 mm in a body-centred cubic (BCC) arrangement.
The porosity is 0.417, close to a real pellet bed (about 0.38 to 0.42).
The BCC lattice has mirror planes, so I only mesh a quarter of one column
of the lattice and put **symmetry planes** on the four sides. This makes the
mesh four times smaller with no loss of information (for steady, symmetric
flow).

**Why shrink the spheres to 95 %.** In a real packing spheres touch at single
points. A mesh cannot resolve a point contact: the cells around it become
infinitely thin. The usual fixes in the literature are to shrink the spheres,
enlarge them (overlap), or add small "bridges". I shrank them so there is a
small gap. This slightly lowers the pressure drop near contacts.

**Meshing with snappyHexMesh.**
1. `blockMesh` makes a box of cubes (background mesh).
2. `snappyHexMesh` refines the cubes near each sphere (level 2 = cells split
   twice, so 4 times smaller), removes the cubes inside the spheres
   (castellation), then moves the points onto the sphere surface (snapping).
3. `checkMesh` checks quality (non-orthogonality, skewness). Here: max
   non-orthogonality about 37 degrees, "Mesh OK".
4. Check that I did it right: the fluid volume of the mesh agrees with the
   exact geometry within 0.02 %.

**How I remove the entrance and exit effects.** The gas needs about one
period to "feel" the bed, and the exit also disturbs the flow. So I run two
beds, 2 and 4 periods long, and use only the **difference**: the extra
pressure drop of the 2 extra periods divided by their length. The same trick
works for the temperature.

**Mesh study.** I repeated two Reynolds numbers on a finer mesh (level 3,
cells half the size). If the answer changes little, the mesh is fine enough.
See the README for the numbers.

**Flow regime.** Reynolds numbers 1 to 100, laminar and steady. Real
induration runs at higher Reynolds numbers (a few hundred to about 1000),
where the flow between pellets is unsteady. That would need a transient
solver (`pimpleFoam`) or a turbulence model, and much more computing time.

## 5. Case 2b: heat transfer from gas to pellets

Gas enters at a dimensionless temperature T = 1, the pellet surfaces are kept
at T = 0. With the flow field frozen, `scalarTransportFoam` solves the
convection-diffusion equation for T (diffusivity = nu / Pr). Far from the
inlet the mixed gas temperature falls exponentially:

    T_out(N2) / T_out(N1) = exp( - h a (N2 - N1) s / (rho cp u) )

so I get the heat transfer coefficient h, and Nu = h d / k. I compare it with
the Wakao-Kaguei and Gunn correlations.

**Honest limits I must say myself:**
- Below Re of about 20 the gas reaches the pellet temperature within one
  pellet layer. The outlet temperature is then extremely small and axial heat
  conduction matters, so the simple formula above does not hold. I show the
  heat result only for Re >= 20.
- Wakao-Kaguei and Gunn are for random beds of real particles, measured in
  experiments with a dispersion correction. An ordered BCC bed with fixed
  surface temperature is not the same system, so I expect differences.
- Constant gas properties, no radiation, pellets at a fixed temperature (no
  conduction inside them). In my induration model the pellet temperature
  changes with time and the gas is hot; that coupling is the next step.

## 6. Questions I may get in an interview, and short answers

**Why OpenFOAM and not Fluent?** It is free, so I could build and share
everything. The physics and the workflow (geometry, mesh, boundary
conditions, solver settings, convergence, mesh study, validation) are the
same in Fluent and STAR-CCM+. The porous coefficients map one-to-one
(`1/alpha = d`, `C2 = f`). I have not used Fluent yet and I am ready to learn
it.

**How do you know your simulation is right?** Three checks: (1) the porous
zone reproduces Ergun exactly (verification), (2) the mesh volume equals the
exact geometry, (3) a finer mesh changes the result by only a few percent.
Then I compare with independent data (Ergun, Wakao-Kaguei, Gunn) and explain
the differences.

**Why is the resolved pressure drop lower than Ergun at higher Re?** The BCC
bed is ordered and aligned with the flow, so there are straighter paths than
in a random bed and less form drag. The 95 % shrink also opens the contacts.
At low Re the result is close to Ergun's viscous term.

**What would you do next?** Hot gas with temperature-dependent properties
(`rhoPorousSimpleFoam` or `buoyantSimpleFoam`); a two-temperature porous model
(gas and solid temperatures separately, as in my induration model); a random
packing made with DEM (for example LIGGGHTS or Blender) instead of BCC; higher
Re with a transient solver; and running in parallel on a cluster.

**What did you NOT do?** No turbulence models, no transient runs, no
combustion or radiation, no compressible flow, no commercial CFD code. I know
the theory of most of these from my studies and my own models, but I have
not run them in a CFD code yet.

## 7. Exercises to really learn it (about 1 to 2 hours each)

1. Change `EPS` to 0.36 in `case1_porous_ergun/run_sweep.py` and re-run. Does
   the pressure drop follow the (1 - eps)^2 / eps^3 trend you expect?
2. In Case 1b, change the wall strip porosity to 0.45 and to 0.50. How does
   the extra gas flow at the wall change?
3. In Case 2, run `python3 make_case.py test --periods 2 --level 1` and then
   `blockMesh`, `snappyHexMesh -overwrite`, `checkMesh`. Open the mesh in
   ParaView and look at the cells near the gaps.
4. Run one Case 2 velocity with `--level 1` and compare the pressure drop with
   level 2. This is a mesh study you did yourself.
5. Open `case2_resolved_bed/runs/L2_N4_Re50` in ParaView (`foamToVTK` first),
   make a slice, colour it by U, and add streamlines.
