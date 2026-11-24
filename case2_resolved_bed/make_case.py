#!/usr/bin/env python3
"""Case 2: write a complete OpenFOAM case for gas flow through a small bed of
real (resolved) spheres.

Geometry: spheres in a body-centred cubic (BCC) arrangement. The flow goes
along z. Because the BCC lattice is mirror-symmetric, we only mesh a quarter
of one lattice column, [0, s/2] x [0, s/2], and put symmetry planes on the four
sides. In each lattice period s there are two layers of spheres:
    layer A: sphere centre on the corner (0, 0)
    layer B: sphere centre on the corner (s/2, s/2)
so every sphere is cut into a quarter by the symmetry planes.

The spheres are shrunk to 95 % of the touching size so that the mesh does not
have to resolve point contacts (a common, documented simplification).

Usage:  python3 make_case.py <case_dir> --periods 4 --velocity 0.05 --level 2
"""
import argparse
import math
import shutil
from pathlib import Path

D = 0.012                      # pellet (sphere) diameter [m]
SHRINK = 0.95                  # d = SHRINK * nearest-neighbour distance
S = D / (SHRINK * math.sqrt(3) / 2)   # BCC lattice constant [m]
NU = 1.516e-5                  # air at 20 C [m2/s]
PR = 0.71                      # Prandtl number of air
L_IN = 1.5 * S                 # empty length before the bed
L_OUT = 3.0 * S                # empty length after the bed

HEADER = "FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"


def porosity():
    return 1 - 2 * (math.pi / 6) * D ** 3 / S ** 3


def write(path, cls, obj, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER.format(cls=cls, obj=obj) + body)


def sphere_centres(periods):
    """Centres of all spheres in the bed (quarter-column coordinates)."""
    h = S / 2
    centres = []
    for k in range(2 * periods):
        z = L_IN + S / 4 + k * S / 2
        xy = (0.0, 0.0) if k % 2 == 0 else (h, h)
        centres.append((xy[0], xy[1], z))
    return centres


def make(case, periods, velocity, level, base_cells):
    case = Path(case)
    if case.exists():
        shutil.rmtree(case)
    h = S / 2
    length = L_IN + periods * S + L_OUT
    dx = h / base_cells
    nz = round(length / dx)

    write(case / "system/blockMeshDict", "dictionary", "blockMeshDict", f"""
// Background mesh: a box of cubes. snappyHexMesh later cuts the spheres out of it.
convertToMeters 1;
vertices
(
    (0 0 0) ({h} 0 0) ({h} {h} 0) (0 {h} 0)
    (0 0 {length}) ({h} 0 {length}) ({h} {h} {length}) (0 {h} {length})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({base_cells} {base_cells} {nz}) simpleGrading (1 1 1) );
boundary
(
    inlet  {{ type patch; faces ((0 3 2 1)); }}
    outlet {{ type patch; faces ((4 5 6 7)); }}
    symX0  {{ type symmetryPlane; faces ((0 4 7 3)); }}
    symX1  {{ type symmetryPlane; faces ((1 2 6 5)); }}
    symY0  {{ type symmetryPlane; faces ((0 1 5 4)); }}
    symY1  {{ type symmetryPlane; faces ((3 7 6 2)); }}
);
""")

    centres = sphere_centres(periods)
    geom = "\n".join(
        f"    pellet{i} {{ type searchableSphere; centre ({x} {y} {z}); radius {D / 2}; }}"
        for i, (x, y, z) in enumerate(centres))
    # a point that is surely in the gas: inside the empty inlet section
    loc = (0.31 * h, 0.17 * h, 0.5 * L_IN)
    write(case / "system/snappyHexMeshDict", "dictionary", "snappyHexMeshDict", f"""
castellatedMesh true;
snap            true;
addLayers       false;

geometry
{{
{geom}
}}

castellatedMeshControls
{{
    maxLocalCells       2000000;
    maxGlobalCells      8000000;
    minRefinementCells  10;
    maxLoadUnbalance    0.10;
    nCellsBetweenLevels 2;
    features            ();
    refinementSurfaces
    {{
        "pellet.*"
        {{
            level ({level} {level});
            patchInfo {{ type wall; inGroups (pellets); }}
        }}
    }}
    resolveFeatureAngle 30;
    refinementRegions   {{}}
    locationInMesh      ({loc[0]} {loc[1]} {loc[2]});
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch    3;
    tolerance       2.0;
    nSolveIter      100;
    nRelaxIter      5;
    nFeatureSnapIter 10;
    implicitFeatureSnap false;
    explicitFeatureSnap false;
    multiRegionFeatureSnap false;
}}

addLayersControls
{{
    relativeSizes true;
    layers {{}}
    expansionRatio 1.0;
    finalLayerThickness 0.3;
    minThickness 0.1;
    nGrow 0;
    featureAngle 60;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}}

meshQualityControls
{{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-30;
    minTetQuality 1e-30;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.05;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
}}

mergeTolerance 1e-6;
""")

    write(case / "system/controlDict", "dictionary", "controlDict", """
application     simpleFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         3000;
deltaT          1;
writeControl    timeStep;
writeInterval   3000;
purgeWrite      1;
writeFormat     ascii;
writePrecision  10;
runTimeModifiable true;
""")
    write(case / "system/fvSchemes", "dictionary", "fvSchemes", """
ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,T)      bounded Gauss linearUpwind grad(T);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited corrected 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited corrected 0.5; }
""")
    write(case / "system/fvSolution", "dictionary", "fvSolution", """
solvers
{
    p  { solver GAMG; smoother GaussSeidel; tolerance 1e-9; relTol 0.05; }
    U  { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.1; }
    T  { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent      yes;
    residualControl { p 1e-6; U 1e-7; T 1e-9; }
    nCorrectors     1;
}
relaxationFactors
{
    fields    { p 1; }
    equations { U 0.8; }
}
""")
    write(case / "system/decomposeParDict", "dictionary", "decomposeParDict", """
numberOfSubdomains 4;
method          scotch;
""")
    write(case / "constant/transportProperties", "dictionary", "transportProperties", f"""
transportModel Newtonian;
nu  {NU};          // air at 20 C [m2/s]
DT  {NU / PR};     // thermal diffusivity of air = nu / Pr [m2/s], used by scalarTransportFoam
""")
    write(case / "constant/turbulenceProperties", "dictionary", "turbulenceProperties", """
// Laminar: we only run particle Reynolds numbers below about 100, where the
// flow between the pellets is laminar and steady.
simulationType laminar;
""")
    sym = "\n".join(f"    {n} {{ type symmetryPlane; }}" for n in ("symX0", "symX1", "symY0", "symY1"))
    write(case / "0/U", "volVectorField", "U", f"""
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 {velocity});
boundaryField
{{
    inlet   {{ type fixedValue; value uniform (0 0 {velocity}); }}
    outlet  {{ type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 {velocity}); }}
    pellets {{ type noSlip; }}
{sym}
}}
""")
    write(case / "0/p", "volScalarField", "p", f"""
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{{
    inlet   {{ type fixedFluxPressure; value uniform 0; }}   // = zero gradient here, but it writes the face values
    outlet  {{ type fixedValue; value uniform 0; }}
    pellets {{ type zeroGradient; }}
{sym}
}}
""")
    # Dimensionless temperature: gas enters at T = 1, pellet surfaces held at T = 0.
    write(case / "0/T", "volScalarField", "T", f"""
dimensions [0 0 0 1 0 0 0];
internalField uniform 0;
boundaryField
{{
    inlet   {{ type fixedValue; value uniform 1; }}
    outlet  {{ type inletOutlet; inletValue uniform 0; value uniform 0; }}
    pellets {{ type fixedValue; value uniform 0; }}
{sym}
}}
""")
    (case / "case_info.txt").write_text(
        f"periods {periods}\nvelocity {velocity}\nlevel {level}\nbase_cells {base_cells}\n"
        f"D {D}\nS {S}\nporosity {porosity()}\nL_in {L_IN}\nL_out {L_OUT}\n")
    return case


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--periods", type=int, default=4)
    ap.add_argument("--velocity", type=float, default=0.05)
    ap.add_argument("--level", type=int, default=2)
    ap.add_argument("--base-cells", type=int, default=10)
    a = ap.parse_args()
    make(a.case, a.periods, a.velocity, a.level, a.base_cells)
    print(f"lattice constant s = {S * 1000:.3f} mm, porosity = {porosity():.4f}")
