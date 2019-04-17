
import numpy as np
from mayavi import mlab
from raysect.core import translate, Point3D
from raysect.optical import World
from raysect.primitive import Sphere, Mesh
from raysect_mayavi.primitives import to_mesh
from raysect_mayavi.primitives.triangle import triangle3d_intersects_triangle3d


world = World()
# s1 = Sphere(0.5, transform=translate(-0.25, 0, 0), name='s1')
s1 = Sphere(0.5, transform=translate(-0.6, 0, 0), name='s1')
# s2 = Sphere(0.5, transform=translate(0.25, 0, 0), name='s2')
s2 = Sphere(0.5, transform=translate(0.6, 0, 0), name='s2')


s1_vertices, s1_triangles = to_mesh(s1)
n_s1_vertices = s1_vertices.shape[0]
n_s1_triangles = s1_triangles.shape[0]
s1_mesh = Mesh(vertices=s1_vertices, triangles=s1_triangles, smoothing=False)

s2_vertices, s2_triangles = to_mesh(s2)
n_s2_vertices = s2_vertices.shape[0]
n_s2_triangles = s2_triangles.shape[0]
s2_mesh = Mesh(vertices=s2_vertices, triangles=s2_triangles, smoothing=False)

combined_vertices = np.zeros((n_s1_vertices+n_s2_vertices, 3))
combined_vertices[0:n_s1_vertices] = s1_vertices[:]
combined_vertices[n_s1_vertices:] = s2_vertices[:]

combined_triangles = np.zeros((n_s1_triangles+n_s2_triangles, 3))
combined_triangles[0:n_s1_triangles] = s1_triangles[:]
combined_triangles[n_s1_triangles:] = s2_triangles[:] + n_s1_vertices

combined_mesh = Mesh(vertices=combined_vertices, triangles=combined_triangles, smoothing=False, parent=world)


s1_intersects = np.zeros(n_s1_triangles)
s2_intersects = np.zeros(n_s2_triangles)
for s1_tri_id in [0]:

    v1, v2, v3 = s1_triangles[s1_tri_id]
    u1x, u1y, u1z = s1_vertices[v1]
    u2x, u2y, u2z = s1_vertices[v2]
    u3x, u3y, u3z = s1_vertices[v3]

    for s2_tri_id in [107]:

        v1, v2, v3 = s2_triangles[s2_tri_id]
        v1x, v1y, v1z = s2_vertices[v1]
        v2x, v2y, v2z = s2_vertices[v2]
        v3x, v3y, v3z = s2_vertices[v3]

        # print('testing', s1_tri_id, s2_tri_id)
        result = triangle3d_intersects_triangle3d(u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
                                                  v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z)

        if result[0]:
            print(True, s1_tri_id, s2_tri_id)
            s1_intersects[s1_tri_id] = 1
            s2_intersects[s2_tri_id] = 1
            line_origin = result[3]
            line_vector = result[4]


dx = s1_vertices[:, 0]
dy = s1_vertices[:, 1]
dz = s1_vertices[:, 2]
# mlab.triangular_mesh(dx, dy, dz, s1_triangles, scalars=s1_intersects)
s = mlab.pipeline.triangular_mesh_source(dx, dy, dz, s1_triangles)
s.data.cell_data.scalars = s1_intersects
surf = mlab.pipeline.surface(s)
surf.contour.filled_contours = True

dx = s2_vertices[:, 0]
dy = s2_vertices[:, 1]
dz = s2_vertices[:, 2]
# mlab.triangular_mesh(dx, dy, dz, s2_triangles, scalars=s2_intersects)
# s = mlab.pipeline.triangular_mesh_source(dx, dy, dz, s2_triangles)
# s.data.cell_data.scalars = s2_intersects
# surf = mlab.pipeline.surface(s)
# surf.contour.filled_contours = True


start = line_origin - line_vector.normalise() * 10
end = line_origin + line_vector.normalise() * 10
print(start, end)
mlab.plot3d([start.x, end.x], [start.y, end.y], [start.z, end.z], tube_radius=0.005)



