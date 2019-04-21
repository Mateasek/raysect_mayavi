
import numpy as np
from mayavi import mlab
from raysect.core import Point3D
from raysect_mayavi.primitives.triangle import Triangle, triangle3d_intersects_triangle3d


# u1 = Point3D(0.0, 0.0, 1.0)
# u2 = Point3D(0.0, 0.0, 0.0)
# u3 = Point3D(1.0, 0.0, 0.0)
#
# v1 = Point3D(0.6000000238418579, 1.600000023841858, 1.600000023841858)
# v2 = Point3D(0.6000000238418579, 0.6000000238418579, 1.600000023841858)
# v3 = Point3D(1.600000023841858, 0.6000000238418579, 1.600000023841858)

# u1 = Point3D(0.0, 0.0, 1.0)
# u2 = Point3D(0.0, 0.0, 0.0)
# u3 = Point3D(1.0, 0.0, 0.0)
#
# v1 = Point3D(0.6000000238418579, 1.600000023841858, 0.6000000238418579)
# v2 = Point3D(0.6000000238418579, 0.6000000238418579, 0.6000000238418579)
# v3 = Point3D(0.6000000238418579, 1.600000023841858, 1.600000023841858)

# u1 = Point3D(0.0, 1.0, 1.0)
# u2 = Point3D(0.0, 0.0, 1.0)
# u3 = Point3D(1.0, 0.0, 1.0)
#
# v1 = Point3D(0.6000000238418579, 0.6000000238418579, 1.600000023841858)
# v2 = Point3D(0.6000000238418579, 0.6000000238418579, 0.6000000238418579)
# v3 = Point3D(1.600000023841858, 0.6000000238418579, 0.6000000238418579)

# u1 = Point3D(0.0, 1.0, 0.0)
# u2 = Point3D(0.0, 0.0, 0.0)
# u3 = Point3D(0.0, 1.0, 1.0)
#
# v1 = Point3D(0.6000000238418579, 0.6000000238418579, 1.600000023841858)
# v2 = Point3D(0.6000000238418579, 0.6000000238418579, 0.6000000238418579)
# v3 = Point3D(1.600000023841858, 0.6000000238418579, 0.6000000238418579)

# u1 = Point3D(0.0, 1.0, 0.0)
# u2 = Point3D(0.0, 0.0, 0.0)
# u3 = Point3D(0.0, 1.0, 1.0)
#
# v1 = Point3D(0.6000000238418579, 1.600000023841858, 1.600000023841858)
# v2 = Point3D(0.6000000238418579, 0.6000000238418579, 1.600000023841858)
# v3 = Point3D(1.600000023841858, 0.6000000238418579, 1.600000023841858)

u1 = Point3D(0.2628655433654785, 0.4253253936767578, 0.0)
u2 = Point3D(0.21694427728652954, 0.4313342273235321, 0.12994596362113953)
u3 = Point3D(0.34689024090766907, 0.35102322697639465, 0.08031101524829865)

v1 = Point3D(0.30000001192092896, -0.5, 0.30000001192092896)
v2 = Point3D(0.30000001192092896, 0.5, -0.699999988079071)
v3 = Point3D(0.30000001192092896, 0.5, 0.30000001192092896)

t1 = Triangle(u1.x, u1.y, u1.z, u2.x, u2.y, u2.z, u3.x, u3.y, u3.z)
t2 = Triangle(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z)

result = triangle3d_intersects_triangle3d(t1, t2)

dx = np.array([u1.x, u2.x, u3.x])
dy = np.array([u1.y, u2.y, u3.y])
dz = np.array([u1.z, u2.z, u3.z])
t1_triangles = np.array([[0, 1, 2]])
mlab.triangular_mesh(dx, dy, dz, t1_triangles)

t1_cp = Point3D(dx.mean(), dy.mean(), dz.mean())
u1u2 = u1.vector_to(u2)
u1u3 = u1.vector_to(u3)
n_pi1 = u1u2.cross(u1u3).normalise()


dx = np.array([v1.x, v2.x, v3.x])
dy = np.array([v1.y, v2.y, v3.y])
dz = np.array([v1.z, v2.z, v3.z])
t2_triangles = np.array([[0, 1, 2]])
mlab.triangular_mesh(dx, dy, dz, t2_triangles)

t2_cp = Point3D(dx.mean(), dy.mean(), dz.mean())
v1v2 = v1.vector_to(v2)
v1v3 = v1.vector_to(v3)
n_pi2 = v1v2.cross(v1v3).normalise()

# Plot plane-1 normal
end = t1_cp + n_pi1.normalise() * 0.5
mlab.plot3d([t1_cp.x, end.x], [t1_cp.y, end.y], [t1_cp.z, end.z], tube_radius=0.005)

# Plot plane-2 normal
end = t2_cp + n_pi2.normalise() * 0.5
mlab.plot3d([t2_cp.x, end.x], [t2_cp.y, end.y], [t2_cp.z, end.z], tube_radius=0.005)


if result[0]:
    t2_point, t3_point = result[1], result[2]
    mlab.plot3d([t2_point.x, t3_point.x], [t2_point.y, t3_point.y], [t2_point.z, t3_point.z], tube_radius=0.005, color=(0, 0, 0))

