

from libc.math cimport abs
from raysect.core.math cimport Vector3D, new_vector3d, Point3D, new_point3d


cpdef tuple triangle3d_intersects_triangle3d(double u1x, double u1y, double u1z,
                                             double u2x, double u2y, double u2z,
                                             double u3x, double u3y, double u3z,
                                             double v1x, double v1y, double v1z,
                                             double v2x, double v2y, double v2z,
                                             double v3x, double v3y, double v3z,
                                             double tolerance=1e-6):
    """
    Cython utility for finding the intersection of two 3d triangles.
    
    If not intersection is found the tuple (False, None) will be returned. If an intersection
    is found, the tuple returned is (True, Intersect).
    
    http://web.stanford.edu/class/cs277/resources/papers/Moller1997b.pdf

    :param double u1x: x coord of triangle u vertex 1.
    :param double u1y: y coord of triangle u vertex 1.
    :param double u1z: z coord of triangle u vertex 1.
    :param double u2x: x coord of triangle u vertex 2.
    :param double u2y: y coord of triangle u vertex 2.
    :param double u2z: z coord of triangle u vertex 2.
    :param double u3x: x coord of triangle u vertex 3.
    :param double u3y: y coord of triangle u vertex 3.
    :param double u3z: z coord of triangle u vertex 3.
    :param double v1x: x coord of triangle v vertex 1.
    :param double v1y: y coord of triangle v vertex 1.
    :param double v1z: z coord of triangle v vertex 1.
    :param double v2x: x coord of triangle v vertex 2.
    :param double v2y: y coord of triangle v vertex 2.
    :param double v2z: z coord of triangle v vertex 2.
    :param double v3x: x coord of triangle v vertex 3.
    :param double v3y: y coord of triangle v vertex 3.
    :param double v3z: z coord of triangle v vertex 3.
    :param double tolerance: the tolerance level of the intersection calculation.
    :return: (False, None) if not intersection is found.
      If an intersection is found, the tuple returned is (True, Intersect).
    :rtype: tuple
    """

    cdef:

        Vector3D u1, u2, u3, v1, v2, v3
        Point3D uc, vc
        Vector3D n_pi1, n_pi2, line_vec
        double d_pi1, d_pi2
        double uv_distance, d_u_from_v, d_v_from_u
        double lo_x, lo_y, lo_z
        double du1, du2, du3, dv1, dv2, dv3
        double pu1, pu2, pu3, pv1, pv2, pv3
        double min_pu, max_pu, min_mv, max_mv
        double t1, t2, t3, t4
        Point3D t2_point, t3_point
        bint valid_intersection_found = False
        bint t1_found = False, t2_found = False
        bint t3_found = False, t4_found = False

    u1 = new_vector3d(u1x, u1y, u1z)
    u2 = new_vector3d(u2x, u2y, u2z)
    u3 = new_vector3d(u3x, u3y, u3z)
    uc = new_point3d((u1x+u2x+u3x)/3, (u1y+u2y+u3y)/3, (u1z+u2z+u3z)/3)

    v1 = new_vector3d(v1x, v1y, v1z)
    v2 = new_vector3d(v2x, v2y, v2z)
    v3 = new_vector3d(v3x, v3y, v3z)
    vc = new_point3d((v1x+v2x+v3x)/3, (v1y+v2y+v3y)/3, (v1z+v2z+v3z)/3)

    # TODO - consider moving to Hessian form of plane equation

    # Setup Plane 1 equation
    n_pi1 = u2.sub(u1).cross(u3.sub(u1))  # normal vector of plane 1
    d_pi1 = n_pi1.dot(u1)  # point satisfying the plane equation for plane 1

    # Setup Plane 2 equation
    n_pi2 = v2.sub(v1).cross(v3.sub(v1))  # normal vector of plane 2
    d_pi2 = n_pi2.dot(v1)  # point satisfying the plane equation for plane 2

    # Find line on intersection of planes P1 and P2
    line_vec = n_pi1.cross(n_pi2)

    # Planes might be parallel
    if line_vec.length == 0:
        if n_pi1.x == n_pi2.x and n_pi1.y == n_pi2.y and n_pi1.z == n_pi2.z and d_pi1 == d_pi2:
            raise NotImplementedError("Planes are parallel and overlapping, this case is not yet implemented.")
        else:
            return False,
    line_vec = line_vec.normalise()

    if n_pi1.y != 0.0:

        denominator = n_pi2.x - (n_pi1.x * n_pi2.y / n_pi1.y)
        if denominator != 0:
            lo_x = (d_pi2 - (d_pi1 * n_pi2.y / n_pi1.y)) / denominator
            lo_y = (d_pi1 - n_pi1.x * lo_x) / n_pi1.y
            lo_z = 0
            valid_intersection_found = True

    if not valid_intersection_found and n_pi2.y != 0.0:

        denominator = n_pi1.x - (n_pi2.x * n_pi1.y / n_pi2.y)
        if denominator != 0:
            lo_x = (d_pi1 - (d_pi2 * n_pi1.y / n_pi2.y)) / denominator
            lo_y = (d_pi2 - n_pi2.x * lo_x) / n_pi2.y
            lo_z = 0
            valid_intersection_found = True

    if not valid_intersection_found and n_pi1.z != 0.0:

        denominator = n_pi2.x - (n_pi1.x * n_pi2.z / n_pi1.z)
        if denominator != 0:
            lo_x = (d_pi2 - (d_pi1 * n_pi2.z / n_pi1.z)) / denominator
            lo_y = 0
            lo_z = (d_pi1 - n_pi1.x * lo_x) / n_pi1.z
            valid_intersection_found = True

    if not valid_intersection_found and n_pi2.z != 0.0:

        denominator = n_pi1.x - (n_pi2.x * n_pi1.z / n_pi2.z)
        if denominator != 0:
            lo_x = (d_pi1 - (d_pi2 * n_pi1.z / n_pi2.z)) / denominator
            lo_y = 0
            lo_z = (d_pi2 - n_pi2.x * lo_x) / n_pi2.z
            valid_intersection_found = True

    if not valid_intersection_found and n_pi1.y != 0.0:

        denominator = n_pi2.z - (n_pi1.z * n_pi2.y / n_pi1.y)
        if denominator != 0:
            lo_x = 0
            lo_z = (d_pi2 - (d_pi1 * n_pi2.y / n_pi1.y)) / denominator
            lo_y = (d_pi1 - n_pi1.z * lo_z) / n_pi1.y
            valid_intersection_found = True

    if not valid_intersection_found and n_pi2.y != 0.0:

        denominator = n_pi1.z - (n_pi2.z * n_pi1.y / n_pi2.y)
        if denominator != 0:
            lo_x = 0
            lo_z = (d_pi1 - (d_pi2 * n_pi1.y / n_pi2.y)) / denominator
            lo_y = (d_pi2 - n_pi2.z * lo_z) / n_pi2.y
            valid_intersection_found = True

    if not valid_intersection_found:
        print('')
        print('Debugging information')
        print('')
        print('u1 = Point3D({}, {}, {})'.format(u1x, u1y, u1z))
        print('u2 = Point3D({}, {}, {})'.format(u2x, u2y, u2z))
        print('u3 = Point3D({}, {}, {})'.format(u3x, u3y, u3z))
        print('')
        print('v1 = Point3D({}, {}, {})'.format(v1x, v1y, v1z))
        print('v2 = Point3D({}, {}, {})'.format(v2x, v2y, v2z))
        print('v3 = Point3D({}, {}, {})'.format(v3x, v3y, v3z))
        raise ValueError("Unsolvable triangle intersection problem.")

    line_origin = new_vector3d(lo_x, lo_y, lo_z)

    l = n_pi2.length
    du1 = (n_pi2.dot(u1) - d_pi2) / l
    du2 = (n_pi2.dot(u2) - d_pi2) / l
    du3 = (n_pi2.dot(u3) - d_pi2) / l

    # case for no intersection
    if (du1 > 0 and du2 > 0 and du3 > 0) or (du1 < 0 and du2 < 0 and du3 < 0):
        return False,

    l = n_pi1.length
    dv1 = (n_pi1.dot(v1) - d_pi1) / l
    dv2 = (n_pi1.dot(v2) - d_pi1) / l
    dv3 = (n_pi1.dot(v3) - d_pi1) / l

    if (dv1 > 0 and dv2 > 0 and dv3 > 0) or (dv1 < 0 and dv2 < 0 and dv3 < 0):
        return False,

    # case for co-planar triangles
    elif (du1 == 0 and du2 == 0 and du3 == 0) or (dv1 == 0 and dv2 == 0 and dv3 == 0):
        raise NotImplementedError("Planes are parallel and overlapping, this case is not yet implemented.")

    # case for overlapping 3D triangles
    else:

        pu1 = line_vec.dot(u1.sub(line_origin))
        pu2 = line_vec.dot(u2.sub(line_origin))
        pu3 = line_vec.dot(u3.sub(line_origin))

        min_pu = min(pu1, pu2, pu3)
        max_pu = max(pu1, pu2, pu3)

        pv1 = line_vec.dot(v1.sub(line_origin))
        pv2 = line_vec.dot(v2.sub(line_origin))
        pv3 = line_vec.dot(v3.sub(line_origin))

        min_pv = min(pv1, pv2, pv3)
        max_pv = max(pv1, pv2, pv3)

        if not du1 - du2 == 0:
            t1 = pu1 + (pu2 - pu1) * (du1 / (du1 - du2))
            if min_pu <= t1 <= max_pu:
                t1_found = True

        if not du1 - du3 == 0:
            if not t1_found:
                t1 = pu1 + (pu3 - pu1) * (du1 / (du1 - du3))
                if min_pu <= t1 <= max_pu:
                    t1_found = True
            else:
                t2 = pu1 + (pu3 - pu1) * (du1 / (du1 - du3))
                if min_pu <= t2 <= max_pu:
                    t2_found = True

        if not t2_found and not du2 - du3 == 0:
            t2 = pu2 + (pu3 - pu2) * (du2 / (du2 - du3))
            if min_pu <= t2 <= max_pu:
                t2_found = True

        # ignore case of single contact point, need two contact points for valid intersection
        if not (t1_found and t2_found):
            return False,

        if t1 > t2:
            t1, t2 = t2, t1

        if not dv1 - dv2 == 0:
            t3 = pv1 + (pv2 - pv1) * (dv1 / (dv1 - dv2))
            if min_pv <= t3 <= max_pv:
                t3_found = True

        if not dv1 - dv3 == 0:
            if not t3_found:
                t3 = pv1 + (pv3 - pv1) * (dv1 / (dv1 - dv3))
                if min_pv <= t3 <= max_pv:
                    t3_found = True
            else:
                t4 = pv1 + (pv3 - pv1) * (dv1 / (dv1 - dv3))
                if min_pv <= t4 <= max_pv:
                    t4_found = True

        if not t4_found and not dv2 - dv3 == 0:
            t4 = pv2 + (pv3 - pv2) * (dv2 / (dv2 - dv3))
            if min_pv <= t4 <= max_pv:
                t4_found = True

        # ignore case of single contact point, need two contact points for valid intersection
        if not (t3_found and t4_found):
            return False,

        if t3 > t4:
            t3, t4 = t4, t3

        # ensure triangles are ordered lowest to highers in terms of parameter t (i.e. left to right)
        if t3 < t1:
            t1, t3 = t3, t1
            t2, t4 = t4, t2

        # test for no intersection
        if (t1 < t3 and t1 < t4 and t2 < t3 and t2 < t4) or (t1 > t3 and t1 > t4 and t2 > t3 and t2 > t4):
            return False,
        else:

            t2_point = new_point3d(line_origin.x + t2 * line_vec.x, line_origin.y + t2 * line_vec.y, line_origin.z + t2 * line_vec.z)
            t3_point = new_point3d(line_origin.x + t3 * line_vec.x, line_origin.y + t3 * line_vec.y, line_origin.z + t3 * line_vec.z)

            if t2_point.distance_to(t3_point) < tolerance:
                return False,

            return True, t2_point, t3_point
