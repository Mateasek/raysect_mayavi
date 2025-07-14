
import numpy as np
from mayavi import mlab
from raysect.core import translate, Point3D
from raysect.optical import World
from raysect.primitive import Sphere, Mesh
from raycanvas.primitives import to_mesh
from raycanvas.primitives.triangle import triangle3d_intersects_triangle3d


world = World()
s1 = Sphere(0.5, transform=translate(-0.25, 0, 0), name='s1')
s2 = Sphere(0.5, transform=translate(0.25, 0, 0), name='s2')


s1_vertices, s1_triangles = to_mesh(s1)
n_s1_vertices = s1_vertices.shape[0]
n_s1_triangles = s1_triangles.shape[0]
s1_mesh = Mesh(vertices=s1_vertices, triangles=s1_triangles, smoothing=False)

s2_vertices, s2_triangles = to_mesh(s2)
n_s2_vertices = s2_vertices.shape[0]
n_s2_triangles = s2_triangles.shape[0]
s2_mesh = Mesh(vertices=s2_vertices, triangles=s2_triangles, smoothing=False)


dx = s1_vertices[:, 0]
dy = s1_vertices[:, 1]
dz = s1_vertices[:, 2]
# mlab.triangular_mesh(dx, dy, dz, s1_triangles, scalars=s1_intersects)
s = mlab.pipeline.triangular_mesh_source(dx, dy, dz, s1_triangles)
s1_highlight = np.zeros(n_s1_triangles)
# s1_highlight[296] = 1
s.data.cell_data.scalars = s1_highlight
surf = mlab.pipeline.surface(s)
surf.contour.filled_contours = True

dx = s2_vertices[:, 0]
dy = s2_vertices[:, 1]
dz = s2_vertices[:, 2]
# mlab.triangular_mesh(dx, dy, dz, s2_triangles, scalars=s2_intersects)
s = mlab.pipeline.triangular_mesh_source(dx, dy, dz, s2_triangles)
s2_highlight = np.zeros(n_s2_triangles)
s2_highlight[216] = 1
s.data.cell_data.scalars = s2_highlight
surf = mlab.pipeline.surface(s)
surf.contour.filled_contours = True

