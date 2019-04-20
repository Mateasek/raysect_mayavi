
import numpy as np
from scipy.spatial import Delaunay, KDTree

from raysect.core cimport Point3D, translate, rotate_basis

from raysect_mayavi.primitives.triangle cimport triangle3d_intersects_triangle3d


cdef class CSG_Operator:

    cdef bint _m1(self, double signed_distance):
        raise NotImplementedError

    cdef bint _m2(self, double signed_distance):
        raise NotImplementedError


cdef class Union(CSG_Operator):

    cdef bint _m1(self, double signed_distance):
        if signed_distance >= 0:
            return True
        return False

    cdef bint _m2(self, double signed_distance):
        if signed_distance >= 0:
            return True
        return False


cdef class Intersect(CSG_Operator):

    cdef bint _m1(self, double signed_distance):
        if signed_distance <= 0:
            return True
        return False

    cdef bint _m2(self, double signed_distance):
        if signed_distance <= 0:
            return True
        return False


cdef class Subtract(CSG_Operator):

    cdef bint _m1(self, double signed_distance):
        if signed_distance >= 0:
            return True
        return False

    cdef bint _m2(self, double signed_distance):
        if signed_distance <= 0:
            return True
        return False


cpdef Mesh perform_mesh_csg(Mesh mesh_1, Mesh mesh_2, CSG_Operator operator):

    cdef:
        tuple result

    # extract vertex and triangle data for mesh 1
    m1_vertices = mesh_1.data.vertices.copy()
    m1_triangles = mesh_1.data.triangles.copy()
    n_m1_vertices = m1_vertices.shape[0]
    n_m1_triangles = m1_triangles.shape[0]

    # calculate triangle centres and KDtree for mesh 1
    m1_tri_centres = np.zeros(m1_triangles.shape)
    for m1_tri_id in range(n_m1_triangles):
        v1, v2, v3 = m1_triangles[m1_tri_id]
        u1x, u1y, u1z = m1_vertices[v1]
        u2x, u2y, u2z = m1_vertices[v2]
        u3x, u3y, u3z = m1_vertices[v3]
        m1_tri_centres[m1_tri_id, :] = (u1x + u2x + u3x) / 3, (u1y + u2y + u3y) / 3, (u1z + u2z + u3z) / 3
    m1_kdtree = KDTree(m1_tri_centres)

    # extract vertex and triangle data for mesh 2
    m2_vertices = mesh_2.data.vertices.copy()
    m2_triangles = mesh_2.data.triangles.copy()
    n_m2_vertices = m2_vertices.shape[0]
    n_m2_triangles = m2_triangles.shape[0]

    # calculate triangle centres and KDtree for mesh 2
    m2_tri_centres = np.zeros(m2_triangles.shape)
    for m2_tri_id in range(n_m2_triangles):
        v1, v2, v3 = m2_triangles[m2_tri_id]
        v1x, v1y, v1z = m2_vertices[v1]
        v2x, v2y, v2z = m2_vertices[v2]
        v3x, v3y, v3z = m2_vertices[v3]
        m2_tri_centres[m2_tri_id, :] = (v1x + v2x + v3x) / 3, (v1y + v2y + v3y) / 3, (v1z + v2z + v3z) / 3
    m2_kdtree = KDTree(m2_tri_centres)

    ####################################################################################################################
    # FIND ALL INTERSECTIONS AND SIGNED DISTANCES

    m1_intersects = np.full((n_m1_triangles), False)
    m1_intersections = {}
    m1_from_m2_signed_distance = np.empty(n_m1_triangles)
    m1_from_m2_signed_distance[:] = 1E999

    m2_intersects = np.full((n_m2_triangles), False)
    m2_intersections = {}
    m2m1_distance = np.empty(n_m2_triangles)
    m2m1_distance[:] = 1E999
    m2_from_m1_signed_distance = np.empty(n_m2_triangles)
    m2_from_m1_signed_distance[:] = 1E999

    for m1_tri_id in range(n_m1_triangles):

        # unpack vertices for triangle in mesh 1
        v1, v2, v3 = m1_triangles[m1_tri_id]
        u1x, u1y, u1z = m1_vertices[v1]
        u1 = Point3D(u1x, u1y, u1z)
        u2x, u2y, u2z = m1_vertices[v2]
        u2 = Point3D(u2x, u2y, u2z)
        u3x, u3y, u3z = m1_vertices[v3]
        u3 = Point3D(u3x, u3y, u3z)
        uc = Point3D(m1_tri_centres[m1_tri_id, 0], m1_tri_centres[m1_tri_id, 1], m1_tri_centres[m1_tri_id, 2])
        n_pi1 = u1.vector_to(u2).cross(u1.vector_to(u3)).normalise()  # normal vector of plane 1
        tri1_radius = max(uc.distance_to(u1), uc.distance_to(u2), uc.distance_to(u3))

        # find all mesh 2 triangles within radius r of mesh 1
        m2_candidates = m2_kdtree.query_ball_point([uc.x, uc.y, uc.z], tri1_radius)

        m1m2_distance = 1E999  # temp variable
        for m2_tri_id in m2_candidates:

            # unpack vertices for triangle in mesh 2
            v1, v2, v3 = m2_triangles[m2_tri_id]
            v1x, v1y, v1z = m2_vertices[v1]
            v1 = Point3D(v1x, v1y, v1z)
            v2x, v2y, v2z = m2_vertices[v2]
            v2 = Point3D(v2x, v2y, v2z)
            v3x, v3y, v3z = m2_vertices[v3]
            v3 = Point3D(v3x, v3y, v3z)
            vc = Point3D(m2_tri_centres[m2_tri_id, 0], m2_tri_centres[m2_tri_id, 1], m2_tri_centres[m2_tri_id, 2])
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

                intersect_pt_1 = result[1]
                intersect_pt_2 = result[2]

                m1_intersects[m1_tri_id] = True
                try:
                    m1_intersections[m1_tri_id].append((m2_tri_id, intersect_pt_1, intersect_pt_2))
                except KeyError:
                    m1_intersections[m1_tri_id] = [(m2_tri_id, intersect_pt_1, intersect_pt_2)]

                m2_intersects[m2_tri_id] = True
                try:
                    m2_intersections[m2_tri_id].append((m1_tri_id, intersect_pt_1, intersect_pt_2))
                except KeyError:
                    m2_intersections[m2_tri_id] = [(m1_tri_id, intersect_pt_1, intersect_pt_2)]

            else:
                if uv_distance < m1m2_distance:
                    m1m2_distance = uv_distance
                    m1_from_m2_signed_distance[m1_tri_id] = d_u_from_v

        if len(m2_candidates) == 0:

            # query closest triangle in mesh 2
            closest_m2_tri = m2_kdtree.query([uc.x, uc.y, uc.z])[1]

            # unpack vertices for closest triangle in mesh 2
            v1, v2, v3 = m2_triangles[closest_m2_tri]
            v1x, v1y, v1z = m2_vertices[v1]
            v1 = Point3D(v1x, v1y, v1z)
            v2x, v2y, v2z = m2_vertices[v2]
            v2 = Point3D(v2x, v2y, v2z)
            v3x, v3y, v3z = m2_vertices[v3]
            v3 = Point3D(v3x, v3y, v3z)
            vc = Point3D(m2_tri_centres[closest_m2_tri, 0], m2_tri_centres[closest_m2_tri, 1], m2_tri_centres[closest_m2_tri, 2])
            n_pi2 = v1.vector_to(v2).cross(v1.vector_to(v3)).normalise()  # normal vector of plane 2

            d_u_from_v = -n_pi2.dot(uc.vector_to(vc))
            m1_from_m2_signed_distance[m1_tri_id] = d_u_from_v


    for m2_tri_id in range(n_m2_triangles):

        # unpack vertices for triangle in mesh 2
        v1, v2, v3 = m2_triangles[m2_tri_id]
        v1x, v1y, v1z = m2_vertices[v1]
        v1 = Point3D(v1x, v1y, v1z)
        v2x, v2y, v2z = m2_vertices[v2]
        v2 = Point3D(v2x, v2y, v2z)
        v3x, v3y, v3z = m2_vertices[v3]
        v3 = Point3D(v3x, v3y, v3z)
        vc = Point3D(m2_tri_centres[m2_tri_id, 0], m2_tri_centres[m2_tri_id, 1], m2_tri_centres[m2_tri_id, 2])
        n_pi2 = v1.vector_to(v2).cross(v1.vector_to(v3)).normalise()  # normal vector of plane 2

        # query closest triangle in mesh 1
        closest_m1_tri = m1_kdtree.query([vc.x, vc.y, vc.z])[1]

        # unpack vertices for closest triangle in mesh 1
        v1, v2, v3 = m1_triangles[closest_m1_tri]
        u1x, u1y, u1z = m1_vertices[v1]
        u1 = Point3D(u1x, u1y, u1z)
        u2x, u2y, u2z = m1_vertices[v2]
        u2 = Point3D(u2x, u2y, u2z)
        u3x, u3y, u3z = m1_vertices[v3]
        u3 = Point3D(u3x, u3y, u3z)
        uc = Point3D(m1_tri_centres[closest_m1_tri, 0], m1_tri_centres[closest_m1_tri, 1], m1_tri_centres[closest_m1_tri, 2])
        n_pi1 = u1.vector_to(u2).cross(u1.vector_to(u3)).normalise()  # normal vector of plane 1

        d_v_from_u = -n_pi1.dot(vc.vector_to(uc))
        m2_from_m1_signed_distance[m2_tri_id] = d_v_from_u


    ####################################################################################################################
    # ADD TRIANGLES THAT DON"T NEED TO BE SPLIT


    combined_vertices = m1_vertices.tolist()
    combined_vertices.extend(m2_vertices.tolist())

    combined_triangles = []
    for m1_tri_id in range(n_m1_triangles):
        if not m1_intersects[m1_tri_id] and operator._m1(m1_from_m2_signed_distance[m1_tri_id]):
            combined_triangles.append(m1_triangles[m1_tri_id, :].tolist())
    for m2_tri_id in range(n_m2_triangles):
        if not m2_intersects[m2_tri_id] and operator._m2(m2_from_m1_signed_distance[m2_tri_id]):
            combined_triangles.append((m2_triangles[m2_tri_id, :] + n_m1_vertices).tolist())


    ####################################################################################################################
    # SPLIT INTERSECTING TRIANGLES


    for m1_tri_id in m1_intersections:

        # unpack triangle 1

        v1, v2, v3 = m1_triangles[m1_tri_id]
        u1x, u1y, u1z = m1_vertices[v1]
        u1 = Point3D(u1x, u1y, u1z)
        u2x, u2y, u2z = m1_vertices[v2]
        u2 = Point3D(u2x, u2y, u2z)
        u3x, u3y, u3z = m1_vertices[v3]
        u3 = Point3D(u3x, u3y, u3z)
        uc = Point3D((u1x + u2x + u3x) / 3, (u1y + u2y + u3y) / 3, (u1z + u2z + u3z) / 3)
        u1u2 = u1.vector_to(u2).normalise()
        u1u3 = u1.vector_to(u3).normalise()
        n_pi1 = u1u2.cross(u1u3).normalise()  # normal vector of plane 1

        # transform to x-y plane
        m1_transform = translate(u1.x, u1.y, u1.z) * rotate_basis(u1u2, n_pi1)
        m1_inv_transform = m1_transform.inverse()

        points = set()
        for vertex in (u1, u2, u3):
            tv = vertex.transform(m1_inv_transform)
            points.add((tv.x, tv.z))

        debug_intersection_points = []
        for intersection in m1_intersections[m1_tri_id]:

            m2_tri_id, ut1, ut2 = intersection

            debug_intersection_points.append((ut1, ut2))

            ut1t = ut1.transform(m1_inv_transform)
            points.add((ut1t.x, ut1t.z))
            ut2t = ut2.transform(m1_inv_transform)
            points.add((ut2t.x, ut2t.z))

            # TODO - replace this with getting the closest triangle in original mesh
            v1, v2, v3 = m2_triangles[m2_tri_id]
            v1x, v1y, v1z = m2_vertices[v1]
            v1 = Point3D(v1x, v1y, v1z)
            v2x, v2y, v2z = m2_vertices[v2]
            v2 = Point3D(v2x, v2y, v2z)
            v3x, v3y, v3z = m2_vertices[v3]
            v3 = Point3D(v3x, v3y, v3z)
            vc = Point3D((v1x + v2x + v3x) / 3, (v1y + v2y + v3y) / 3, (v1z + v2z + v3z) / 3)
            n_pi2 = v1.vector_to(v2).cross(v1.vector_to(v3)).normalise()  # normal vector of plane 2

        # transform unique points back to world space
        world_points = []
        points = list(points)
        num_current_vertices = len(combined_vertices)
        for vertex in points:
            vertex = Point3D(vertex[0], 0, vertex[1])
            tv = vertex.transform(m1_transform)
            combined_vertices.append([tv.x, tv.y, tv.z])

        # perform Delaunay triangulation
        vertices_2d = np.array(points)
        new_triangles = Delaunay(vertices_2d).simplices

        debug_new_triangles = []
        for new_tri_id in range(new_triangles.shape[0]):

            # unpack triangle 1
            v1, v2, v3 = new_triangles[new_tri_id] + num_current_vertices
            u1x, u1y, u1z = combined_vertices[v1]
            u2x, u2y, u2z = combined_vertices[v2]
            u3x, u3y, u3z = combined_vertices[v3]
            uc = Point3D((u1x + u2x + u3x) / 3, (u1y + u2y + u3y) / 3, (u1z + u2z + u3z) / 3)

            d_u_from_v = -n_pi2.dot(uc.vector_to(vc))

            debug_new_triangles.append([[u1x, u1y, u1z], [u2x, u2y, u2z], [u3x, u3y, u3z], [u1x, u1y, u1z]])

            if operator._m1(d_u_from_v):
                combined_triangles.append([v1, v2, v3])


    for m2_tri_id in m2_intersections:

        # unpack triangle 2

        v1, v2, v3 = m2_triangles[m2_tri_id]
        v1x, v1y, v1z = m2_vertices[v1]
        v1 = Point3D(v1x, v1y, v1z)
        v2x, v2y, v2z = m2_vertices[v2]
        v2 = Point3D(v2x, v2y, v2z)
        v3x, v3y, v3z = m2_vertices[v3]
        v3 = Point3D(v3x, v3y, v3z)
        vc = Point3D((v1x + v2x + v3x) / 3, (v1y + v2y + v3y) / 3, (v1z + v2z + v3z) / 3)
        v1v2 = v1.vector_to(v2).normalise()
        v1v3 = v1.vector_to(v3).normalise()
        n_pi2 = v1v2.cross(v1v3).normalise()  # normal vector of plane 2

        # transform to x-y plane
        m2_transform = translate(v1.x, v1.y, v1.z) * rotate_basis(v1v2, n_pi2)
        m2_inv_transform = m2_transform.inverse()

        points = set()
        for vertex in (v1, v2, v3):
            tv = vertex.transform(m2_inv_transform)
            points.add((tv.x, tv.z))

        debug_intersection_points = []
        for intersection in m2_intersections[m2_tri_id]:

            m1_tri_id, ut1, ut2 = intersection

            ut1t = ut1.transform(m2_inv_transform)
            points.add((ut1t.x, ut1t.z))
            ut2t = ut2.transform(m2_inv_transform)
            points.add((ut2t.x, ut2t.z))

            # TODO - replace this with getting the closest triangle in original mesh
            v1, v2, v3 = m1_triangles[m1_tri_id]
            u1x, u1y, u1z = m1_vertices[v1]
            u1 = Point3D(u1x, u1y, u1z)
            u2x, u2y, u2z = m1_vertices[v2]
            u2 = Point3D(u2x, u2y, u2z)
            u3x, u3y, u3z = m1_vertices[v3]
            u3 = Point3D(u3x, u3y, u3z)
            uc = Point3D((u1x + u2x + u3x) / 3, (u1y + u2y + u3y) / 3, (u1z + u2z + u3z) / 3)
            u1u2 = u1.vector_to(u2).normalise()
            u1u3 = u1.vector_to(u3).normalise()
            n_pi1 = u1u2.cross(u1u3).normalise()  # normal vector of plane 1

        # transform unique points back to world space
        world_points = []
        points = list(points)
        num_current_vertices = len(combined_vertices)
        for vertex in points:
            vertex = Point3D(vertex[0], 0, vertex[1])
            tv = vertex.transform(m2_transform)
            combined_vertices.append([tv.x, tv.y, tv.z])

        # perform Delaunay triangulation
        vertices_2d = np.array(points)
        new_triangles = Delaunay(vertices_2d).simplices

        for new_tri_id in range(new_triangles.shape[0]):

            # unpack triangle 1
            v1, v2, v3 = new_triangles[new_tri_id] + num_current_vertices
            v1x, v1y, v1z = combined_vertices[v1]
            v2x, v2y, v2z = combined_vertices[v2]
            v3x, v3y, v3z = combined_vertices[v3]
            vc = Point3D((v1x + v2x + v3x) / 3, (v1y + v2y + v3y) / 3, (v1z + v2z + v3z) / 3)

            d_v_from_u = -n_pi1.dot(vc.vector_to(uc))

            if operator._m2(d_v_from_u):
                combined_triangles.append([v1, v2, v3])

    return Mesh(combined_vertices, combined_triangles)

