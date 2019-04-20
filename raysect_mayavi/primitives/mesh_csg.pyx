
import numpy as np
from libc.math cimport sqrt
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


cpdef Mesh perform_mesh_csg(Mesh mesh_1, Mesh mesh_2, CSG_Operator operator, double tolerance=1e-6):

    cdef:
        tuple result

    # extract vertex and triangle data for mesh 1
    m1_vertices = mesh_1.data.vertices.copy()
    m1_triangles = mesh_1.data.triangles.copy()
    n_m1_vertices = m1_vertices.shape[0]
    n_m1_triangles = m1_triangles.shape[0]

    # extract vertex and triangle data for mesh 2
    m2_vertices = mesh_2.data.vertices.copy()
    m2_triangles = mesh_2.data.triangles.copy()
    n_m2_vertices = m2_vertices.shape[0]
    n_m2_triangles = m2_triangles.shape[0]

    # calculate triangle centres, longest side length and KDtree for mesh 2
    m2_tri_centres = np.zeros(m2_triangles.shape)
    longest_m2_side = 0
    for m2_tri_id in range(n_m2_triangles):
        v1, v2, v3 = m2_triangles[m2_tri_id]
        v1x, v1y, v1z = m2_vertices[v1]
        v2x, v2y, v2z = m2_vertices[v2]
        v3x, v3y, v3z = m2_vertices[v3]
        m2_tri_centres[m2_tri_id, :] = (v1x + v2x + v3x) / 3, (v1y + v2y + v3y) / 3, (v1z + v2z + v3z) / 3
        v1v2_length = sqrt((v2x-v1x)**2+(v2y-v1y)**2+(v2z-v1z)**2)
        v2v3_length = sqrt((v3x-v2x)**2+(v3y-v2y)**2+(v3z-v2z)**2)
        v1v3_length = sqrt((v3x-v1x)**2+(v3y-v1y)**2+(v3z-v1z)**2)
        longest_m2_side = max(longest_m2_side, v1v2_length, v2v3_length, v1v3_length)
    m2_kdtree = KDTree(m2_tri_centres)

    ####################################################################################################################
    # FIND ALL INTERSECTIONS

    m1_intersects = np.full((n_m1_triangles), False)
    m1_split_points = {}

    m2_intersects = np.full((n_m2_triangles), False)
    m2_split_points = {}

    for m1_tri_id in range(n_m1_triangles):

        # unpack vertices for triangle in mesh 1
        v1, v2, v3 = m1_triangles[m1_tri_id]
        u1x, u1y, u1z = m1_vertices[v1]
        u1 = Point3D(u1x, u1y, u1z)
        u2x, u2y, u2z = m1_vertices[v2]
        u2 = Point3D(u2x, u2y, u2z)
        u3x, u3y, u3z = m1_vertices[v3]
        u3 = Point3D(u3x, u3y, u3z)
        uc = Point3D((u1x + u2x + u3x) / 3, (u1y + u2y + u3y) / 3, (u1z + u2z + u3z) / 3)
        n_pi1 = u1.vector_to(u2).cross(u1.vector_to(u3)).normalise()  # normal vector of plane 1

        # find all mesh 2 triangles within radius r of mesh 1
        m2_candidates = m2_kdtree.query_ball_point([uc.x, uc.y, uc.z], longest_m2_side)

        for m2_tri_id in m2_candidates:

            # unpack vertices for triangle in mesh 2
            v1, v2, v3 = m2_triangles[m2_tri_id]
            v1x, v1y, v1z = m2_vertices[v1]
            v1 = Point3D(v1x, v1y, v1z)
            v2x, v2y, v2z = m2_vertices[v2]
            v2 = Point3D(v2x, v2y, v2z)
            v3x, v3y, v3z = m2_vertices[v3]
            v3 = Point3D(v3x, v3y, v3z)
            vc = Point3D((v1x + v2x + v3x) / 3, (v1y + v2y + v3y) / 3, (v1z + v2z + v3z) / 3)
            n_pi2 = v1.vector_to(v2).cross(v1.vector_to(v3)).normalise()  # normal vector of plane 2

            uv_distance = uc.distance_to(vc)
            d_u_from_v = -n_pi2.dot(uc.vector_to(vc))
            d_v_from_u = -n_pi1.dot(vc.vector_to(uc))

            try:
                result = triangle3d_intersects_triangle3d(u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
                                                          v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
                                                          tolerance=tolerance)
            except NotImplementedError:
                continue

            if result[0]:

                intsct_pt_1 = result[1]
                intsct_pt_2 = result[2]

                m1_intersects[m1_tri_id] = True
                try:
                    m1_split_points[m1_tri_id]
                except KeyError:
                    m1_split_points[m1_tri_id] = [u1, u2, u3]

                # add intersection points if separation from other split points is greater than tolerance
                distances = [intsct_pt_1.distance_to(pt) for pt in m1_split_points[m1_tri_id]]
                if all([d > tolerance for d in distances]):
                    m1_split_points[m1_tri_id].append(intsct_pt_1)
                distances = [intsct_pt_2.distance_to(pt) for pt in m1_split_points[m1_tri_id]]
                if all([d > tolerance for d in distances]):
                    m1_split_points[m1_tri_id].append(intsct_pt_2)

                m2_intersects[m2_tri_id] = True
                try:
                    m2_split_points[m2_tri_id]
                except KeyError:
                    m2_split_points[m2_tri_id] = [v1, v2, v3]

                # add intersection points if separation from other split points is greater than tolerance
                distances = [intsct_pt_1.distance_to(pt) for pt in m2_split_points[m2_tri_id]]
                if all([d > tolerance for d in distances]):
                    m2_split_points[m2_tri_id].append(intsct_pt_1)
                distances = [intsct_pt_2.distance_to(pt) for pt in m2_split_points[m2_tri_id]]
                if all([d > tolerance for d in distances]):
                        m2_split_points[m2_tri_id].append(intsct_pt_2)


    ####################################################################################################################
    # MAKE NEW MESHES

    # extract vertex and triangle data for mesh 1
    m1_new_vertices = m1_vertices.tolist()
    m1_new_triangles = []
    for m1_tri_id in range(n_m1_triangles):
        if not m1_intersects[m1_tri_id]:
            m1_new_triangles.append(m1_triangles[m1_tri_id, :].tolist())

    m2_new_vertices = m2_vertices.tolist()
    m2_new_triangles = []
    for m2_tri_id in range(n_m2_triangles):
        if not m2_intersects[m2_tri_id]:
            m2_new_triangles.append(m2_triangles[m2_tri_id, :].tolist())


    ####################################################################################################################
    # SPLIT INTERSECTING TRIANGLES AND ADD THEM TO THE NEW MESH

    for m1_tri_id in m1_split_points:

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

        num_current_m1_vertices = len(m1_new_vertices)
        vertices_2d = []
        for point in m1_split_points[m1_tri_id]:
            tv = point.transform(m1_inv_transform)
            vertices_2d.append([tv.x, tv.z])
            m1_new_vertices.append([point.x, point.y, point.z])

        # perform Delaunay triangulation
        vertices_2d = np.array(vertices_2d)
        new_triangles = Delaunay(vertices_2d).simplices

        for new_tri_id in range(new_triangles.shape[0]):

            v1, v2, v3 = new_triangles[new_tri_id]
            u1 = m1_split_points[m1_tri_id][v1]
            u2 = m1_split_points[m1_tri_id][v2]
            u3 = m1_split_points[m1_tri_id][v3]
            u1u2 = u1.vector_to(u2).normalise()
            u1u3 = u1.vector_to(u3).normalise()
            n_new = u1u2.cross(u1u3).normalise()  # normal vector of plane 1

            if not n_pi1.dot(n_new) > 0:
                new_triangles[new_tri_id, 1], new_triangles[new_tri_id, 2] = new_triangles[new_tri_id, 2], new_triangles[new_tri_id, 1]

        m1_new_triangles.extend(new_triangles + num_current_m1_vertices)


    for m2_tri_id in m2_split_points:

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

        num_current_m2_vertices = len(m2_new_vertices)
        vertices_2d = []
        for point in m2_split_points[m2_tri_id]:
            tv = point.transform(m2_inv_transform)
            vertices_2d.append([tv.x, tv.z])
            m2_new_vertices.append([point.x, point.y, point.z])

        # perform Delaunay triangulation
        vertices_2d = np.array(vertices_2d)
        new_triangles = Delaunay(vertices_2d).simplices

        for new_tri_id in range(new_triangles.shape[0]):

            v1, v2, v3 = new_triangles[new_tri_id]
            u1 = m2_split_points[m2_tri_id][v1]
            u2 = m2_split_points[m2_tri_id][v2]
            u3 = m2_split_points[m2_tri_id][v3]
            u1u2 = u1.vector_to(u2).normalise()
            u1u3 = u1.vector_to(u3).normalise()
            n_new = u1u2.cross(u1u3).normalise()  # normal vector of plane 1

            if not n_pi2.dot(n_new) > 0:
                new_triangles[new_tri_id, 1], new_triangles[new_tri_id, 2] = new_triangles[new_tri_id, 2], new_triangles[new_tri_id, 1]

        m2_new_triangles.extend(new_triangles + num_current_m2_vertices)

    ####################################################################################################################
    # Recompute mesh KDTrees

    # extract vertex and triangle data for mesh 1
    m1_vertices = np.array(m1_new_vertices)
    m1_triangles = np.array(m1_new_triangles)
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
    m2_vertices = np.array(m2_new_vertices)
    m2_triangles = np.array(m2_new_triangles)
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
    # Perform Signed distance operations

    combined_vertices = m1_vertices.tolist()
    combined_vertices.extend(m2_vertices.tolist())

    combined_triangles = []

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

        if operator._m1(d_u_from_v):
            combined_triangles.append(m1_triangles[m1_tri_id, :].tolist())

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

        if operator._m2(d_v_from_u):
            combined_triangles.append((m2_triangles[m2_tri_id, :] + n_m1_vertices).tolist())

    return Mesh(combined_vertices, combined_triangles)

