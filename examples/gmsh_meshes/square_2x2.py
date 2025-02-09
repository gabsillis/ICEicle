import gmsh
import math
import sys
import numpy as np


# mesh size parameter
lc = 0.05

# setup the gmsh model
gmsh.initialize()
gmsh.model.add("square_2x2")

# =========================
# = make the bounding box =
# =========================
bp1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
bp2 = gmsh.model.geo.addPoint(2, 0, 0, lc)
bp3 = gmsh.model.geo.addPoint(2, 2, 0, lc)
bp4 = gmsh.model.geo.addPoint(0, 2, 0, lc)


bl1 = gmsh.model.geo.addLine(bp1, bp2);
bl2 = gmsh.model.geo.addLine(bp2, bp3);
bl3 = gmsh.model.geo.addLine(bp3, bp4);
bl4 = gmsh.model.geo.addLine(bp4, bp1);
bbox_loop = gmsh.model.geo.addCurveLoop([bl1, bl2,  bl3, bl4])
domain = gmsh.model.geo.addPlaneSurface([bbox_loop])

# =======================
# = Boundary Conditions =
# =======================
gmsh.model.addPhysicalGroup(1, [bl1], 1) # initial condition
gmsh.model.addPhysicalGroup(1, [bl2], 2) # right
gmsh.model.addPhysicalGroup(1, [bl3], 3) # spacetime future
gmsh.model.addPhysicalGroup(1, [bl4], 4) # left 
gmsh.model.addPhysicalGroup(2, [domain], 5)

gmsh.model.geo.synchronize()
gmsh.option.setNumber("Mesh.Algorithm", 5)  # delunay
gmsh.model.mesh.generate(2)
gmsh.write("square_2x2.msh");

gmsh.fltk.run()

gmsh.finalize()
