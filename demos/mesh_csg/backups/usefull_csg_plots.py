
import numpy as np
from scipy.spatial import Delaunay
from mayavi import mlab

from raysect.core import translate, Point3D, rotate_basis
from raysect.optical import World
from raysect.primitive import Sphere, Mesh
from raysect_mayavi.primitives import to_mesh
from raysect_mayavi.primitives.triangle import triangle3d_intersects_triangle3d


class Union:

    def t1(self, x):
        if x >= 0:
            return True
        return False

    def t2(self, x):
        if x >= 0:
            return True
        return False


class Intersection:

    def t1(self, x):
        if x <= 0:
            return True
        return False

    def t2(self, x):
        if x <= 0:
            return True
        return False


class Difference:

    def t1(self, x):
        if x >= 0:
            return True
        return False

    def t2(self, x):
        if x <= 0:
            return True
        return False


# OPERATION = Union()
OPERATION = Intersection()
# OPERATION = Difference()


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


########################################################################################################################
# FIND ALL INTERSECTIONS AND SIGNED DISTANCES

s1_intersects = np.full((n_s1_triangles), False)
s1_intersections = {}
s1s2_distance = np.empty(n_s1_triangles)
s1s2_distance[:] = 1E999
s1_from_s2_signed_distance = np.empty(n_s1_triangles)
s1_from_s2_signed_distance[:] = 1E999

s2_intersects = np.full((n_s2_triangles), False)
s2_intersections = {}
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

        try:
            result = triangle3d_intersects_triangle3d(u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
                                                      v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z)
        except NotImplementedError:
            continue

        if result[0]:

            ut1 = result[1]
            ut2 = result[2]

            s1_intersects[s1_tri_id] = True
            try:
                s1_intersections[s1_tri_id].append((s2_tri_id, ut1, ut2))
            except KeyError:
                s1_intersections[s1_tri_id] = [(s2_tri_id, ut1, ut2)]

            s2_intersects[s2_tri_id] = True
            try:
                s2_intersections[s2_tri_id].append((s1_tri_id, ut1, ut2))
            except KeyError:
                s2_intersections[s2_tri_id] = [(s1_tri_id, ut1, ut2)]

        else:
            if uv_distance < s1s2_distance[s1_tri_id]:
                s1s2_distance[s1_tri_id] = uv_distance
                s1_from_s2_signed_distance[s1_tri_id] = d_u_from_v

            if uv_distance < s2s1_distance[s2_tri_id]:
                s2s1_distance[s2_tri_id] = uv_distance
                s2_from_s1_signed_distance[s2_tri_id] = d_v_from_u


########################################################################################################################
# ADD TRIANGLES THAT DON"T NEED TO BE SPLIT


combined_vertices = s1_vertices.tolist()
combined_vertices.extend(s2_vertices.tolist())

combined_triangles = []
for s1_tri_id in range(n_s1_triangles):
    if not s1_intersects[s1_tri_id] and OPERATION.t1(s1_from_s2_signed_distance[s1_tri_id]):
        combined_triangles.append(s1_triangles[s1_tri_id, :].tolist())
for s2_tri_id in range(n_s2_triangles):
    if not s2_intersects[s2_tri_id] and OPERATION.t2(s2_from_s1_signed_distance[s2_tri_id]):
        combined_triangles.append((s2_triangles[s2_tri_id, :] + n_s1_vertices).tolist())


########################################################################################################################
# SPLIT INTERSECTING TRIANGLES


for s1_tri_id in s1_intersections:

    print()
    print('s1_tri_id', s1_tri_id, 'len(s1_intersections[s1_tri_id])', len(s1_intersections[s1_tri_id]))

    # unpack triangle 1

    v1, v2, v3 = s1_triangles[s1_tri_id]
    u1x, u1y, u1z = s1_vertices[v1]
    u1 = Point3D(u1x, u1y, u1z)
    u2x, u2y, u2z = s1_vertices[v2]
    u2 = Point3D(u2x, u2y, u2z)
    u3x, u3y, u3z = s1_vertices[v3]
    u3 = Point3D(u3x, u3y, u3z)
    uc = Point3D((u1x + u2x + u3x) / 3, (u1y + u2y + u3y) / 3, (u1z + u2z + u3z) / 3)
    u1u2 = u1.vector_to(u2).normalise()
    u1u3 = u1.vector_to(u3).normalise()
    n_pi1 = u1u2.cross(u1u3).normalise()  # normal vector of plane 1

    debug_vertices = [[u1x, u1y, u1z], [u2x, u2y, u2z], [u3x, u3y, u3z]]
    debug_triangles = [[0, 1, 2]]

    # transform to x-y plane
    s1_transform = translate(u1.x, u1.y, u1.z) * rotate_basis(u1u2, n_pi1)
    s1_inv_transform = s1_transform.inverse()

    points = set()
    for vertex in (u1, u2, u3):
        tv = vertex.transform(s1_inv_transform)
        points.add((tv.x, tv.z))

    debug_intersection_points = []
    for intersection in s1_intersections[s1_tri_id]:

        s2_tri_id, ut1, ut2 = intersection
        print('-->')
        print('s2_tri_id', s2_tri_id)
        print(ut1.distance_to(ut2))

        debug_intersection_points.append((ut1, ut2))

        ut1t = ut1.transform(s1_inv_transform)
        points.add((ut1t.x, ut1t.z))
        ut2t = ut2.transform(s1_inv_transform)
        points.add((ut2t.x, ut2t.z))

        # TODO - replace this with getting the closest triangle in original mesh
        v1, v2, v3 = s2_triangles[s2_tri_id]
        v1x, v1y, v1z = s2_vertices[v1]
        v1 = Point3D(v1x, v1y, v1z)
        v2x, v2y, v2z = s2_vertices[v2]
        v2 = Point3D(v2x, v2y, v2z)
        v3x, v3y, v3z = s2_vertices[v3]
        v3 = Point3D(v3x, v3y, v3z)
        vc = Point3D((v1x + v2x + v3x) / 3, (v1y + v2y + v3y) / 3, (v1z + v2z + v3z) / 3)
        n_pi2 = v1.vector_to(v2).cross(v1.vector_to(v3)).normalise()  # normal vector of plane 2

        offset = len(debug_vertices)
        debug_vertices.extend([[v1x, v1y, v1z], [v2x, v2y, v2z], [v3x, v3y, v3z]])
        debug_triangles.extend([[0+offset, 1+offset, 2+offset]])

    # transform unique points back to world space
    world_points = []
    points = list(points)
    num_current_vertices = len(combined_vertices)
    for vertex in points:
        vertex = Point3D(vertex[0], 0, vertex[1])
        tv = vertex.transform(s1_transform)
        combined_vertices.append([tv.x, tv.y, tv.z])

    # perform Delaunay triangulation
    vertices_2d = np.array(points)
    new_triangles = Delaunay(vertices_2d).simplices

    debug_new_triangles = []
    for new_tri_id in range(new_triangles.shape[0]):

        # unpack triangle 1
        v1, v2, v3 = new_triangles[new_tri_id] + num_current_vertices
        u1x, u1y, u1z = combined_vertices[v1]
        u1 = Point3D(u1x, u1y, u1z)
        u2x, u2y, u2z = combined_vertices[v2]
        u2 = Point3D(u2x, u2y, u2z)
        u3x, u3y, u3z = combined_vertices[v3]
        u3 = Point3D(u3x, u3y, u3z)
        uc = Point3D((u1x + u2x + u3x) / 3, (u1y + u2y + u3y) / 3, (u1z + u2z + u3z) / 3)

        d_u_from_v = -n_pi2.dot(uc.vector_to(vc))

        debug_new_triangles.append([[u1x, u1y, u1z], [u2x, u2y, u2z], [u3x, u3y, u3z], [u1x, u1y, u1z]])

        if OPERATION.t1(d_u_from_v):
            combined_triangles.append([v1, v2, v3])

    # debug_vertices = np.array(debug_vertices)
    # debug_triangles = np.array(debug_triangles)
    #
    # dx = debug_vertices[:, 0]
    # dy = debug_vertices[:, 1]
    # dz = debug_vertices[:, 2]
    # mlab.figure(1)
    # mlab.clf()
    # mlab.triangular_mesh(dx, dy, dz, debug_triangles)
    # for pair in debug_intersection_points:
    #     mlab.plot3d([pair[0].x, pair[1].x], [pair[0].y, pair[1].y], [pair[0].z, pair[1].z], tube_radius=0.001, color=(0, 0, 1))
    #
    # mlab.figure(2)
    # mlab.clf()
    # mlab.triangular_mesh(dx, dy, dz, debug_triangles)
    # for new_tri in debug_new_triangles:
    #     new_tri = np.array(new_tri)
    #     mlab.plot3d(new_tri[:, 0], new_tri[:, 1], new_tri[:, 2], tube_radius=0.005, color=(0, 0, 1))
    #
    # input('...')



combined_vertices = np.array(combined_vertices)
combined_triangles = np.array(combined_triangles)

dx = combined_vertices[:, 0]
dy = combined_vertices[:, 1]
dz = combined_vertices[:, 2]
mlab.triangular_mesh(dx, dy, dz, combined_triangles)


