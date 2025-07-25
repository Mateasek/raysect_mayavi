from raycanvas import visualise_scenegraph

from raysect.core import translate
from raysect.optical import World
from raysect.primitive import Cylinder
from raysect.primitive.lens.spherical import BiConvex, BiConcave, PlanoConvex, PlanoConcave, Meniscus

from mayavi import mlab
# Display lens and cylinder primitives

world = World()
diameter = 1
center_thickness = 0.5
front_curvature = 2
back_curvature = 2

biconvex_lens = BiConvex(diameter, center_thickness, front_curvature, back_curvature, parent=world)
biconcave_lens = BiConcave(diameter, center_thickness, front_curvature, back_curvature, parent=world, transform=translate(0, 0, 1))
planoconvex_lens = PlanoConvex(diameter, center_thickness, front_curvature, parent=world, transform=translate(0, 0, 2))
planoconcave_lens = PlanoConcave(diameter, center_thickness, front_curvature, parent=world, transform=translate(0, 0, 3))
meniscus_lens = Meniscus(diameter, center_thickness, front_curvature, back_curvature, parent=world, transform=translate(0, 0, 4))
cylinder_primitive = Cylinder(radius=0.5 * diameter, height=center_thickness, parent=world, transform=translate(0, 0, 5))


visualise_scenegraph(world)
mlab.show()