
import numpy as np
from mayavi import mlab
from raysect.core import translate, Point3D
from raysect.optical import World
from raysect.primitive import Sphere, Mesh
from raysect_mayavi.primitives import to_mesh
from raysect_mayavi.primitives.triangle import triangle3d_intersects_triangle3d


# OPERATION = 'UNION'
# OPERATION = 'INTERSECTION'
OPERATION = 'DIFFERENCE'


world = World()
s1 = Sphere(0.5, transform=translate(-0.25, 0, 0), name='s1')
s2 = Sphere(0.5, transform=translate(0.25, 0, 0), name='s2')


s1_vertices, s1_triangles = to_mesh(s1)
n_s1_vertices = s1_vertices.shape[0]
n_s1_triangles = s1_triangles.shape[0]
s1_mesh = Mesh(vertices=s1_vertices, triangles=s1_triangles, smoothing=False)
print()
print('n_s1_triangles', n_s1_triangles)

s2_vertices, s2_triangles = to_mesh(s2)
n_s2_vertices = s2_vertices.shape[0]
n_s2_triangles = s2_triangles.shape[0]
s2_mesh = Mesh(vertices=s2_vertices, triangles=s2_triangles, smoothing=False)
print('n_s2_triangles', n_s2_triangles)


s1_intersects = np.zeros(n_s1_triangles)
s1s2_distance = np.empty(n_s1_triangles)
s1s2_distance[:] = 1E999
s1_from_s2_signed_distance = np.empty(n_s1_triangles)
s1_from_s2_signed_distance[:] = 1E999

s2_intersects = np.zeros(n_s2_triangles)
s2s1_distance = np.empty(n_s2_triangles)
s2s1_distance[:] = 1E999
s2_from_s1_signed_distance = np.empty(n_s2_triangles)
s2_from_s1_signed_distance[:] = 1E999
s2_pair = np.zeros(n_s1_triangles)


for s1_tri_id in range(n_s1_triangles):

    v1, v2, v3 = s1_triangles[s1_tri_id]
    u1x, u1y, u1z = s1_vertices[v1]
    u1 = Point3D(u1x, u1y, u1z)
    u2x, u2y, u2z = s1_vertices[v2]
    u2 = Point3D(u2x, u2y, u2z)
    u3x, u3y, u3z = s1_vertices[v3]
    u3 = Point3D(u3x, u3y, u3z)
    uc = Point3D((u1x + u2x + u3x) / 3, (u1y + u2y + u3y) / 3, (u1z + u2z + u3z) / 3)
    n_pi1 = u1.vector_to(u2).cross(u1.vector_to(u3)).normalise()  # normal vector of plane 1

    for s2_tri_id in range(n_s2_triangles):

        v1, v2, v3 = s2_triangles[s2_tri_id]
        v1x, v1y, v1z = s2_vertices[v1]
        v1 = Point3D(v1x, v1y, v1z)
        v2x, v2y, v2z = s2_vertices[v2]
        v2 = Point3D(v2x, v2y, v2z)
        v3x, v3y, v3z = s2_vertices[v3]
        v3 = Point3D(v3x, v3y, v3z)
        vc = Point3D((v1x + v2x + v3x) / 3, (v1y + v2y + v3y) / 3, (v1z + v2z + v3z) / 3)
        n_pi2 = v1.vector_to(v2).cross(v1.vector_to(v3)).normalise()  # normal vector of plane 2

        uv_distance = uc.distance_to(vc)
        d_u_from_v = -n_pi2.dot(uc.vector_to(vc))
        d_v_from_u = -n_pi1.dot(vc.vector_to(uc))

        if uv_distance < s1s2_distance[s1_tri_id]:
            s1s2_distance[s1_tri_id] = uv_distance
            s1_from_s2_signed_distance[s1_tri_id] = d_u_from_v

        if uv_distance < s2s1_distance[s2_tri_id]:
            s2s1_distance[s2_tri_id] = uv_distance
            s2_from_s1_signed_distance[s2_tri_id] = d_v_from_u


if OPERATION == 'UNION':

    combined_vertices = np.zeros((n_s1_vertices+n_s2_vertices, 3))
    combined_vertices[0:n_s1_vertices] = s1_vertices[:]
    combined_vertices[n_s1_vertices:] = s2_vertices[:]

    n_unioned_s1_triangles = sum(d >= 0 for d in s1_from_s2_signed_distance)
    n_unioned_s2_triangles = sum(d >= 0 for d in s2_from_s1_signed_distance)

    combined_triangles = np.zeros((n_unioned_s1_triangles+n_unioned_s2_triangles, 3))

    tri_i = 0
    for s1_tri_id in range(n_s1_triangles):
        if s1_from_s2_signed_distance[s1_tri_id] >= 0:
            combined_triangles[tri_i, :] = s1_triangles[s1_tri_id, :]
            tri_i += 1
    for s2_tri_id in range(n_s2_triangles):
        if s2_from_s1_signed_distance[s2_tri_id] >= 0:
            # print('s2 id', s2_tri_id, 's1 pair', s2_pair[s2_tri_id])
            combined_triangles[tri_i, :] = s2_triangles[s2_tri_id, :] + n_s1_vertices
            tri_i += 1

elif OPERATION == 'INTERSECTION':

    combined_vertices = np.zeros((n_s1_vertices+n_s2_vertices, 3))
    combined_vertices[0:n_s1_vertices] = s1_vertices[:]
    combined_vertices[n_s1_vertices:] = s2_vertices[:]

    n_unioned_s1_triangles = sum(d <= 0 for d in s1_from_s2_signed_distance)
    n_unioned_s2_triangles = sum(d <= 0 for d in s2_from_s1_signed_distance)

    combined_triangles = np.zeros((n_unioned_s1_triangles+n_unioned_s2_triangles, 3))

    tri_i = 0
    for s1_tri_id in range(n_s1_triangles):
        if s1_from_s2_signed_distance[s1_tri_id] <= 0:
            combined_triangles[tri_i, :] = s1_triangles[s1_tri_id, :]
            tri_i += 1
    for s2_tri_id in range(n_s2_triangles):
        if s2_from_s1_signed_distance[s2_tri_id] <= 0:
            # print('s2 id', s2_tri_id, 's1 pair', s2_pair[s2_tri_id])
            combined_triangles[tri_i, :] = s2_triangles[s2_tri_id, :] + n_s1_vertices
            tri_i += 1

elif OPERATION == 'DIFFERENCE':

    combined_vertices = np.zeros((n_s1_vertices+n_s2_vertices, 3))
    combined_vertices[0:n_s1_vertices] = s1_vertices[:]
    combined_vertices[n_s1_vertices:] = s2_vertices[:]

    n_unioned_s1_triangles = sum(d >= 0 for d in s1_from_s2_signed_distance)
    n_unioned_s2_triangles = sum(d <= 0 for d in s2_from_s1_signed_distance)

    combined_triangles = np.zeros((n_unioned_s1_triangles+n_unioned_s2_triangles, 3))

    tri_i = 0
    for s1_tri_id in range(n_s1_triangles):
        if s1_from_s2_signed_distance[s1_tri_id] >= 0:
            combined_triangles[tri_i, :] = s1_triangles[s1_tri_id, :]
            tri_i += 1
    for s2_tri_id in range(n_s2_triangles):
        if s2_from_s1_signed_distance[s2_tri_id] <= 0:
            # print('s2 id', s2_tri_id, 's1 pair', s2_pair[s2_tri_id])
            combined_triangles[tri_i, :] = s2_triangles[s2_tri_id, :] + n_s1_vertices
            tri_i += 1


dx = combined_vertices[:, 0]
dy = combined_vertices[:, 1]
dz = combined_vertices[:, 2]
mlab.triangular_mesh(dx, dy, dz, combined_triangles)


