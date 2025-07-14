import numpy as np

from raysect.primitive import Cone, Mesh

from raycanvas.backend.base.source import TriangularMeshSource
from raycanvas.backend.builtin.primitives.weld_vertices import weld_vertices
from raycanvas.backend.builtin.primitives.geometric.utility import triangulate_cylinder_lid, triangulate_conical_cap


class ConeSource(TriangularMeshSource):
    """
    Class for graphical representation of the Raysect Cone primitive.
    :param raysect_object: Raysect Cone primitive instance
    :param subdivision_count: Number of subdivisions to perform on the basic icosohedren with 12 vertices and 20 faces
    """

    def __init__(self, raysect_object, vertical_divisions=10, cylindrical_divisions=36, radial_divisions=5):

        if not isinstance(raysect_object, Cone):
            raise TypeError("The raysect_object has to be instance of Raysect Cone primitive, wrong type '{}' given.".format(type(raysect_object)))

        self._cap_vertices = None
        self._cap_triangles = None
        self._cone_vertices = None
        self._cone_triangles = None

        self.cylindrical_divisions = cylindrical_divisions
        self.radial_divisions = radial_divisions
        self.vertical_divisions = vertical_divisions

        super().__init__(raysect_object)

    def _graphic_source_from_raysect_object(self):

        if self._cap_vertices is None or self._cap_triangles is None:
            self._generate_cap_surface()
        if self._cone_vertices is None or self._cone_triangles is None:
            self._generate_conic_surface()

        #combine vertices of the 2 surfaces
        vertices = np.concatenate((self._cap_vertices, self._cone_vertices))
        # recalculate and combine triangle ids for merged vertices
        n_cap = self._cap_vertices.shape[0]

        parabola_triangles = self._cone_triangles.copy()
        parabola_triangles += n_cap

        triangles = np.concatenate((self._cap_triangles, parabola_triangles))

        mesh = Mesh(vertices=vertices, triangles=triangles)
        self._raysect_mesh = weld_vertices(mesh)

    def _generate_cap_surface(self):
        """
        Generate vertices and triangles for the cylindrical cap of the parabolic primitive
        :return:
        """

        self._cap_vertices, self._cap_triangles = triangulate_cylinder_lid(self._raysect_object.radius,
                                                                           self.radial_divisions,
                                                                           cylindrical_divisions=self.cylindrical_divisions)

    def _generate_conic_surface(self):

        if self._cap_vertices is None:
            self._generate_cap_surface()

        edge = self._cap_vertices[0:self.cylindrical_divisions, 0:2]
        self._cone_vertices, self._cone_triangles = triangulate_conical_cap(self._raysect_object.height,
                                                                                      self.raysect_object.radius,
                                                                                      self.vertical_divisions,
                                                                                      edge_vertices=edge)