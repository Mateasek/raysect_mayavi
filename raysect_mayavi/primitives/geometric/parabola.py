from raysect_mayavi.primitives.source import SourceMesh
from raysect_mayavi.primitives.geometric.utility import triangulate_open_cylinder, triangulate_cylinder_lid, triangulate_spherical_cap
from raysect_mayavi.primitives.weld_vertices import weld_vertices

from raysect.primitive import Mesh
from raysect.primitive import Cylinder
from raysect.primitive.lens.spherical import BiConvex, BiConcave, PlanoConvex, PlanoConcave, Meniscus

import numpy as np
from math import sqrt


class CylindricalSource(SourceMesh):
    """
    Class for creting mesh sources for objects with cylindrical symmetry e.g. cylinder and lenses. The mesh is composed
    of 3 surfaces: front surface, back surface and an open cylinder surface.
    """

    def __init__(self, raysect_primitive, vertical_divisions=10, cylindrical_divisions=36, radial_divisions=5):

        self._cap_vertices = None
        self._cap_triangles = None
        self._parabola_vertices = None
        self._parabola_triangles = None

        self.cylindrical_divisions = cylindrical_divisions
        self.radial_divisions = radial_divisions

        super().__init__(raysect_primitive)

    def _mesh_from_primitive(self):

        if self._cap_vertices is None or self._cap_triangles is None:
            self._generate_cap_surface()
        if self._parabola_vertices is None or self._parabola_triangles is None:
            self._generate_parabola_surface()

        #combine vertices of the 2 surfaces
        vertices = np.concatenate((self._cap_vertices, self._parabola_vertices))
        # recalculate and combine triangle ids for merged vertices
        n_cap = self._cap_vertices.shape[0]

        parabola_triangles = self._parabola_triangles.copy()
        parabola_triangles += n_cap

        triangles = np.concatenate((self._cap_triangles, parabola_triangles))

        mesh = Mesh(vertices=vertices, triangles=triangles)
        self._raysect_mesh = weld_vertices(mesh)

    def _generate_cap_surface(self):
        """
        Generate vertices and triangles for the cylindrical cap of the parabolic primitive
        :return:
        """

        self._cap_vertices, self._cap_triangles = triangulate_cylinder_lid(self._raysect_primitive.radius,
                                                       self.radial_divisions)

    def _generate_parabola_surface(self):


        radius = self._raysect_primitive.radius
        height = self._raysect_primitive.height


        k = height / (radius * radius)
        a = k * (direction.x * direction.x + direction.y * direction.y)
        b = 2 * k * (direction.x * origin.x + direction.y * origin.y) + direction.z
        c = k * (origin.x * origin.x + origin.y * origin.y) - (height - origin.z)

        edge = self._barrel_vertices[-self.cylindrical_divisions::, 0:2]
