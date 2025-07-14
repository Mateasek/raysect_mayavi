
import numpy as np
from mayavi import mlab
from raysect.core import translate, Point3D, Vector3D
from raysect.optical import World
from raysect.primitive import Sphere, Mesh
from raycanvas.primitives import to_mesh
from raycanvas.primitives.triangle import triangle3d_intersects_triangle3d


def plot_triangles(t1p1, t1p2, t1p3, t2p1, t2p2, t2p3):

    result = triangle3d_intersects_triangle3d(t1p1.x, t1p1.y, t1p1.z, t1p2.x, t1p2.y, t1p2.z, t1p3.x, t1p3.y, t1p3.z,
                                              t2p1.x, t2p1.y, t2p1.z, t2p2.x, t2p2.y, t2p2.z, t2p3.x, t2p3.y, t2p3.z)

    dx = np.array([t1p1.x, t1p2.x, t1p3.x])
    dy = np.array([t1p1.y, t1p2.y, t1p3.y])
    dz = np.array([t1p1.z, t1p2.z, t1p3.z])
    t1_triangles = np.array([[0, 1, 2]])
    mlab.triangular_mesh(dx, dy, dz, t1_triangles)

    t1_cp = Point3D(dx.mean(), dy.mean(), dz.mean())

    dx = np.array([t2p1.x, t2p2.x, t2p3.x])
    dy = np.array([t2p1.y, t2p2.y, t2p3.y])
    dz = np.array([t2p1.z, t2p2.z, t2p3.z])
    t2_triangles = np.array([[0, 1, 2]])
    mlab.triangular_mesh(dx, dy, dz, t2_triangles)

    t2_cp = Point3D(dx.mean(), dy.mean(), dz.mean())

    print('uv_distance', result[1])
    print('d_u_from_v', result[2])
    print('d_v_from_u', result[3])

    if result[0]:

        print('Intersection is True!')

        # Plot plane-1 normal
        n_pi1, d_pi1 = result[8]
        end = t1_cp + n_pi1.normalise() * 0.5
        mlab.plot3d([t1_cp.x, end.x], [t1_cp.y, end.y], [t1_cp.z, end.z], tube_radius=0.005)

        # Plot plane-2 normal
        n_pi2, d_pi2 = result[9]
        end = t2_cp + n_pi2.normalise() * 0.5
        mlab.plot3d([t2_cp.x, end.x], [t2_cp.y, end.y], [t2_cp.z, end.z], tube_radius=0.005)

        # # Plot line along intersection of the two planes
        line_origin = result[6]
        line_vector = result[7]
        # start = line_origin + line_vector * -5
        # end = line_origin + line_vector * -10
        # mlab.plot3d([start.x, end.x], [start.y, end.y], [start.z, end.z], tube_radius=0.005, color=(0, 0, 0))

        # print(line_origin)
        # # solve at z = 0
        # a = np.array([[n_pi1.x, n_pi1.y], [n_pi2.x, n_pi2.y]])
        # b = np.array([d_pi1, d_pi2])
        # x = np.linalg.solve(a, b)
        # print(x)

        # solve at mean triangle 1 x-point
        xs = t1_cp.x
        a = np.array([[n_pi1.y, n_pi1.z], [n_pi2.y, n_pi2.z]])
        b = np.array([d_pi1 - n_pi1.x * xs, d_pi2 - n_pi2.x * xs])
        x = np.linalg.solve(a, b)
        alternative = Point3D(xs, x[0], x[1])
        end = alternative - line_vector
        mlab.plot3d([alternative.x, end.x], [alternative.y, end.y], [alternative.z, end.z], tube_radius=0.005, color=(0, 0, 0))

        # # Plot line along intersection of the two planes
        # line_origin = result[3]
        # line_vector = result[4]
        # end = line_origin + line_vector
        # mlab.plot3d([line_origin.x, end.x], [line_origin.y, end.y], [line_origin.z, end.z], tube_radius=0.005, color=(0, 0, 0))

        # xp = Point3D(cp.x + 0.2, cp.y, (plane1[3] - (plane1[0] * cp.x + plane1[1] * cp.y))/plane1[2])
        # yp = Point3D(cp.x, cp.y + 0.2, (plane1[3] - (plane1[0] * cp.x + plane1[2] * cp.z))/plane1[1])
        # mlab.plot3d([cp.x, xp.x], [cp.y, xp.y], [cp.z, xp.z], tube_radius=0.005)
        # mlab.plot3d([cp.x, yp.x], [cp.y, yp.y], [cp.z, yp.z], tube_radius=0.005)

    else:
        print('Intersection is False!')


# Case 1 - Triangles that do not intersect
# u1 = Point3D(-0.6, 0.2628655560595545, -0.42532540417602765)
# u2 = Point3D(-0.5196889821799917, 0.346890238780218, -0.3510232223880888)
# u3 = Point3D(-0.4700540434961197, 0.21694428227633783, -0.4313342402080971)
#
# v1 = Point3D(0.6, -0.13663326445632884, -0.48096917889196067)
# v2 = Point3D(0.4700540434961197, -0.21694428227633783, -0.4313342402080971)
# v3 = Point3D(0.46856722197021355, -0.08122992405822319, -0.47552825814757654)

# # Case 2 - Triangles that do not intersect
# u1 = Point3D(-0.8628655560595544, 0.42532540417602765, 0.0)
# u2 = Point3D(-0.946890238780218, 0.3510232223880888, 0.08031101782000824)
# u3 = Point3D(-0.8169442822763378, 0.4313342402080971, 0.12994595650388027)
#
# v1 = Point3D(0.8628655560595544, -0.42532540417602765, 0.0)
# v2 = Point3D(0.7366332644563288, -0.48096917889196067, 0.0)
# v3 = Point3D(0.8169442822763378, -0.4313342402080971, -0.12994595650388027)

# # Case 3 - Triangles that do not intersect
# u1 = Point3D(-0.38305571772366215, 0.4313342402080971, 0.12994595650388027)
# u2 = Point3D(-0.4454915028125328, 0.4045084971874727, 0.2500000000000056)
# u3 = Point3D(-0.3061073738537669, 0.3440954801177964, 0.21266270208800983)
#
# v1 = Point3D(0.33713444394042846, 0.0, 0.42532540417601705)
# v2 = Point3D(0.34999999999999437, -0.15450849718746715, 0.4045084971874727)
# v3 = Point3D(0.46856722197021355, -0.08122992405822319, 0.47552825814757654)


# Triangles that do intersect
# u1 = Point3D(0, 0, 0)
# u2 = Point3D(1, 0, 0)
# u3 = Point3D(0, 1, 1)
#
# v1 = Point3D(0, 1, 0) + Vector3D(0.15, 0, 0)
# v2 = Point3D(1, 1, 0) + Vector3D(0.15, 0, 0)
# v3 = Point3D(0, 0, 1) + Vector3D(0.15, 0, 0)

# # Triangles that do intersect
u1 = Point3D(0.012865556059554473, 0.42532540417602765, 0.0)
u2 = Point3D(-0.03305571772366217, 0.4313342402080971, 0.12994595650388027)
u3 = Point3D(0.09689023878021802, 0.3510232223880888, 0.08031101782000824)
v1 = Point3D(-0.012865556059554473, 0.42532540417602765, 0.0)
v2 = Point3D(-0.09689023878021802, 0.3510232223880888, 0.08031101782000824)
v3 = Point3D(0.03305571772366217, 0.4313342402080971, 0.12994595650388027)

# # Triangles just touching, they should not intersect!!!
# Point3D(0.012865556059554473, 0.42532540417602765, 0.0)
# Point3D(-0.03305571772366217, 0.4313342402080971, 0.12994595650388027)
# Point3D(0.09689023878021802, 0.3510232223880888, 0.08031101782000824)
# Point3D(-0.09689023878021802, 0.3510232223880888, 0.08031101782000824)
# Point3D(-0.04389262614623307, 0.3440954801177964, 0.21266270208800983)
# Point3D(0.03305571772366217, 0.4313342402080971, 0.12994595650388027)




result = plot_triangles(u1, u2, u3, v1, v2, v3)

mlab.plot3d([0, 0.5, 1], [0, 0, 0], [0, 0, 0], tube_radius=0.005, color=(1, 0, 0))
mlab.plot3d([0, 0, 0], [0, 0.5, 1], [0, 0, 0], tube_radius=0.005, color=(0, 1, 0))
mlab.plot3d([0, 0, 0], [0, 0, 0], [0, 0.5, 1], tube_radius=0.005, color=(0, 0, 1))


# u1x = -0.6
# u1y = 0.2628655560595545
# u1z = -0.42532540417602765
# u2x = -0.5196889821799917
# u2y = 0.346890238780218
# u2z = -0.3510232223880888
# u3x = -0.4700540434961197
# u3y = 0.21694428227633783
# u3z = -0.4313342402080971
#
# v1x = 0.6
# v1y = -0.13663326445632884
# v1z = -0.48096917889196067
# v2x = 0.4700540434961197
# v2y = -0.21694428227633783
# v2z = -0.4313342402080971
# v3x = 0.46856722197021355
# v3y = -0.08122992405822319
# v3z = -0.47552825814757654
#
#
# result = triangle3d_intersects_triangle3d(u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
#                                           v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z)
#
# line_origin = result[3]
# line_vector = result[4]
# plane1 = result[5]
# plane2 = result[6]
#
# dx = np.array([u1x, u2x, u3x])
# dy = np.array([u1y, u2y, u3y])
# dz = np.array([u1z, u2z, u3z])
# s1_triangles = np.array([[0, 1, 2]])
# # mlab.triangular_mesh(dx, dy, dz, s1_triangles, scalars=s1_intersects)
# s = mlab.pipeline.triangular_mesh_source(dx, dy, dz, s1_triangles)
# s.data.cell_data.scalars = np.array([0])
# surf = mlab.pipeline.surface(s)
# surf.contour.filled_contours = True
#
# cp = Point3D(dx.mean(), dy.mean(), dz.mean())
# xp = Point3D(cp.x + 0.2, cp.y, (plane1[3] - (plane1[0] * cp.x + plane1[1] * cp.y))/plane1[2])
# yp = Point3D(cp.x, cp.y + 0.2, (plane1[3] - (plane1[0] * cp.x + plane1[2] * cp.z))/plane1[1])
# mlab.plot3d([cp.x, xp.x], [cp.y, xp.y], [cp.z, xp.z], tube_radius=0.005)
# mlab.plot3d([cp.x, yp.x], [cp.y, yp.y], [cp.z, yp.z], tube_radius=0.005)
#
#
# dx = np.array([v1x, v2x, v3x])
# dy = np.array([v1y, v2y, v3y])
# dz = np.array([v1z, v2z, v3z])
# s1_triangles = np.array([[0, 1, 2]])
# # mlab.triangular_mesh(dx, dy, dz, s1_triangles, scalars=s1_intersects)
# s = mlab.pipeline.triangular_mesh_source(dx, dy, dz, s1_triangles)
# s.data.cell_data.scalars = np.array([1])
# surf = mlab.pipeline.surface(s)
# surf.contour.filled_contours = True

# start = line_origin + line_vector.normalise() * 5
# end = line_origin + line_vector.normalise() * 15
# mlab.plot3d([start.x, end.x], [start.y, end.y], [start.z, end.z], tube_radius=0.005)


