
import numpy as np
from scipy.spatial import Delaunay
from mayavi import mlab

from raysect.core import translate, Point3D, rotate_basis
from raysect.optical import World
from raysect.primitive import Sphere, Mesh, Intersect, Subtract, Union
from raysect_mayavi.primitives import to_mesh
from raysect_mayavi.primitives.triangle import triangle3d_intersects_triangle3d


world = World()
s1 = Sphere(0.5, transform=translate(-0.25, 0, 0), name='s1')
s2 = Sphere(0.5, transform=translate(0.25, 0, 0), name='s2')
# lens = Intersect(s1, s2, parent=world)
lens = Subtract(s1, s2, parent=world)


from raysect_mayavi import visualise_scenegraph

visualise_scenegraph(world)

