import numpy as np
from mayavi import mlab
from raysect.core import translate, Point3D
from raysect.optical import World
from raysect.primitive import Sphere, Mesh
from raycanvas.backend.builtin.primitives.mesh import to_mesh
from raycanvas.backend.builtin.primitives.triangle import triangle3d_intersects_triangle3d


world = World()
s1 = Sphere(0.5, transform=translate(-0.25, 0, 0), name='s1')
# s1 = Sphere(0.5, transform=translate(-0.6, 0, 0), name='s1')
s2 = Sphere(0.5, transform=translate(0.25, 0, 0), name='s2')
# s2 = Sphere(0.5, transform=translate(0.6, 0, 0), name='s2')


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
# for s1_tri_id in range(n_s1_triangles):
for s1_tri_id in [5]:

    v1, v2, v3 = s1_triangles[s1_tri_id]
    u1x, u1y, u1z = s1_vertices[v1]
    u2x, u2y, u2z = s1_vertices[v2]
    u3x, u3y, u3z = s1_vertices[v3]

    # for s2_tri_id in range(n_s2_triangles):
    for s2_tri_id in [0, 1]:

        v1, v2, v3 = s2_triangles[s2_tri_id]
        v1x, v1y, v1z = s2_vertices[v1]
        v2x, v2y, v2z = s2_vertices[v2]
        v3x, v3y, v3z = s2_vertices[v3]

        try:
            # print('testing', s1_tri_id, s2_tri_id)
            result = triangle3d_intersects_triangle3d(u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z,
                                                      v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z)

            dx = np.array([u1x, u2x, u3x])
            dy = np.array([u1y, u2y, u3y])
            dz = np.array([u1z, u2z, u3z])
            t1_triangles = np.array([[0, 1, 2]])
            mlab.triangular_mesh(dx, dy, dz, t1_triangles)

            dx = np.array([v1x, v2x, v3x])
            dy = np.array([v1y, v2y, v3y])
            dz = np.array([v1z, v2z, v3z])
            t2_triangles = np.array([[0, 1, 2]])
            mlab.triangular_mesh(dx, dy, dz, t2_triangles)

            if result[0]:
                print(True, s1_tri_id, s2_tri_id)
                s1_intersects[s1_tri_id] = 1
                s2_intersects[s2_tri_id] = 1

                ut1 = result[1]
                ut2 = result[2]
                mlab.plot3d([ut1.x, ut2.x], [ut1.y, ut2.y], [ut1.z, ut2.z], tube_radius=0.001,
                            color=(0, 0, 1))

            #     # Plot line along intersection of the two planes
            #     line_origin = result[3]
            #     line_vector = result[4]
            #     end = line_origin + line_vector
            #     mlab.plot3d([line_origin.x, end.x], [line_origin.y, end.y], [line_origin.z, end.z], tube_radius=0.005, color=(0, 0, 0))
            #
            #     print(Point3D(u1x, u1y, u1z))
            #     print(Point3D(u2x, u2y, u2z))
            #     print(Point3D(u3x, u3y, u3z))
            #
            #     print(Point3D(v1x, v1y, v1z))
            #     print(Point3D(v2x, v2y, v2z))
            #     print(Point3D(v3x, v3y, v3z))
            #
            input("waiting...")

        except NotImplementedError:

            continue

            # print('NotImplementedError encountered')
            #
            # dx = np.array([u1x, u2x, u3x])
            # dy = np.array([u1y, u2y, u3y])
            # dz = np.array([u1z, u2z, u3z])
            # t1_triangles = np.array([[0, 1, 2]])
            # mlab.triangular_mesh(dx, dy, dz, t1_triangles)
            #
            # dx = np.array([v1x, v2x, v3x])
            # dy = np.array([v1y, v2y, v3y])
            # dz = np.array([v1z, v2z, v3z])
            # t1_triangles = np.array([[0, 1, 2]])
            # mlab.triangular_mesh(dx, dy, dz, t1_triangles)
            #
            # mlab.plot3d([0, 0.5, 1], [0, 0, 0], [0, 0, 0], tube_radius=0.005, color=(1, 0, 0))
            # mlab.plot3d([0, 0, 0], [0, 0.5, 1], [0, 0, 0], tube_radius=0.005, color=(0, 1, 0))
            # mlab.plot3d([0, 0, 0], [0, 0, 0], [0, 0.5, 1], tube_radius=0.005, color=(0, 0, 1))
            #
            # input("wating...")

# dx = s1_vertices[:, 0]
# dy = s1_vertices[:, 1]
# dz = s1_vertices[:, 2]
# # mlab.triangular_mesh(dx, dy, dz, s1_triangles, scalars=s1_intersects)
# s = mlab.pipeline.triangular_mesh_source(dx, dy, dz, s1_triangles)
# s.data.cell_data.scalars = s1_intersects
# surf = mlab.pipeline.surface(s)
# surf.contour.filled_contours = True


# for vertex_id, vertex in enumerate(s1_vertices):
#     p = Point3D(vertex[0], vertex[1], vertex[2])
# 
#     overlapping_triangles = s2_mesh.data.items_containing(p)
#     if overlapping_triangles:
#         print()
#         print(vertex_id)
#         print(overlapping_triangles)

# from raycanvas import visualise_scenegraph

# visualise_scenegraph(world)

