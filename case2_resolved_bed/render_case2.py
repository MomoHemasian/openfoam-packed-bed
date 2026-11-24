#!/usr/bin/env python3
"""Optional 3D pictures of one resolved run (needs: pip install pyvista).
ParaView gives the same pictures by hand; this script only makes them repeatable."""
import sys
from pathlib import Path

import pyvista as pv

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "plots"
CASE = HERE / "runs" / (sys.argv[1] if len(sys.argv) > 1 else "L2_N4_Re50")


def load():
    (CASE / "case.foam").touch()
    r = pv.POpenFOAMReader(str(CASE / "case.foam"))
    r.set_active_time_value(r.time_values[-1])
    r.enable_all_patch_arrays()
    return r.read()


def mirror4(mesh, h):
    """The run is a quarter of one lattice column. Mirror it in x and y to show the full column."""
    a = pv.merge([mesh, mesh.reflect((1, 0, 0), point=(h, 0, 0))])
    return pv.merge([a, a.reflect((0, 1, 0), point=(0, h, 0))])


def picture(data, field, label, cmap, fname):
    internal = data["internalMesh"]
    h = internal.bounds[1]
    names = [n for n in data["boundary"].keys() if n.startswith("pellet")]
    full = mirror4(internal, h)
    pellets = mirror4(pv.merge([data["boundary"][n] for n in names]), h)
    # Vertical plane y = s/2 through the full column: it cuts the centre spheres
    # through their middle and passes through the open channels between the corner spheres.
    n = (0, 1, 0)
    sl = full.slice(normal=n, origin=(0, h, 0)).cell_data_to_point_data()
    if field == "U":
        sl[label] = (sl.point_data["U"] ** 2).sum(axis=1) ** 0.5
    else:
        sl[label] = sl.point_data[field]
    back = pellets.clip(normal=n, origin=(0, h, 0), invert=False)
    p = pv.Plotter(off_screen=True, window_size=(2400, 640))
    p.set_background("white")
    p.add_mesh(sl, scalars=label, cmap=cmap, lighting=False,
               scalar_bar_args=dict(title=label, color="black", vertical=False, position_x=0.35,
                                    position_y=0.03, width=0.3, height=0.10, fmt="%.2f"))
    p.add_mesh(back, color="#b9b5ad", smooth_shading=True)
    zmid = 0.5 * (internal.bounds[4] + internal.bounds[5])
    fp = (h, h, zmid)
    p.camera.focal_point = fp
    p.camera.position = (h, h - 0.5, zmid)    # look along +y
    p.camera.up = (-1, 0, 0)                    # so the gas flows left to right
    p.camera.parallel_projection = True
    p.camera.parallel_scale = 0.019
    p.add_text("gas flow  -->", position="upper_left", font_size=14, color="black")
    p.screenshot(OUT / fname)


def main():
    data = load()
    picture(data, "U", "|U| (m/s)", "Blues", "case2_velocity.png")
    picture(data, "T", "T (-)", "Oranges", "case2_temperature.png")


if __name__ == "__main__":
    main()
