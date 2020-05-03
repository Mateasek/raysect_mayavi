from raysect.core import Point3D
from raysect.primitive import Mesh

class SourceMesh:

    def __init__(self, raysect_primitive):

        self._raysect_primitive = raysect_primitive
        self._raysect_mesh = None

        self._mesh_from_primitive()

    def _vertices_transform(self):
        vertices = self._raysect_mesh.data.vertices.copy()

        if self._raysect_primitive.parent:
            to_world = self._raysect_primitive.to_root()
        else:
            to_world = self._raysect_primitive.transform

        for i in range(vertices.shape[0]):
            p = Point3D(vertices[i, 0], vertices[i, 1], vertices[i, 2]).transform(to_world)
            vertices[i, 0] = p.x
            vertices[i, 1] = p.y
            vertices[i, 2] = p.z

        return vertices

    @property
    def raysect_primitive(self):
        return self._raysect_primitive

    @property
    def vertices_local(self):
        return self._raysect_mesh.data.vertices()

    @property
    def vertices(self):
        return self._vertices_transform()

    @property
    def vertices_x(self):
        return self._vertices_transform()[:, 0]

    @property
    def vertices_y(self):
        return self._vertices_transform()[:, 1]

    @property
    def vertices_z(self):
        return self._vertices_transform()[:, 2]

    @property
    def triangles(self):
        return self._raysect_mesh.data.triangles

    def _mesh_from_primitive(self):
        """
        Constructs raysect Mesh from sources.
        :return:
        """
        raise NotImplementedError("Virtual method _mesh_from_primitive() has not been implemented.")