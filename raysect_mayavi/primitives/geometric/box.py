from raysect_mayavi.primitives.source import TriangularMeshSource
from raysect_mayavi.primitives.mesh_tools import subdivide

from raysect.core import Point3D
from raysect.primitive import Box, Mesh


class BoxSource(TriangularMeshSource):
    """
    Triangular mesh representation of the Raysect Box primitive.
    :param raysect_object: Raysect Box primitive
    """

    def __init__(self, raysect_object):
        
        if not isinstance(raysect_object, Box):
            raise TypeError("The raysect_object has to be instance of Raysect Box primitive, wrong type '{}' given.".format(type(raysect_object)))
        
        super().__init__(raysect_object)
        
    def _mayavi_source_from_raysect_object(self):
         
        lower = self._raysect_object.lower
        upper = self._raysect_object.upper
        # more negative face in x-z plane
        p1a = lower  # lower corner in x-z plane
        p2a = Point3D(lower.x, lower.y, upper.z)
        p3a = Point3D(upper.x, lower.y, upper.z)  # upper corner in x-z plane
        p4a = Point3D(upper.x, lower.y, lower.z)
        # more positive face in x-z plane
        p1b = Point3D(lower.x, upper.y, lower.z)
        p2b = Point3D(lower.x, upper.y, upper.z)
        p3b = upper
        p4b = Point3D(upper.x, upper.y, lower.z)
        vertices = [[p1a.x, p1a.y, p1a.z], [p2a.x, p2a.y, p2a.z],
                    [p3a.x, p3a.y, p3a.z], [p4a.x, p4a.y, p4a.z],
                    [p1b.x, p1b.y, p1b.z], [p2b.x, p2b.y, p2b.z],
                    [p3b.x, p3b.y, p3b.z], [p4b.x, p4b.y, p4b.z]]
        triangles = [[1, 0, 3], [1, 3, 2],  # front face (x-z)
                    [7, 4, 5], [7, 5, 6],  # rear face (x-z)
                    [5, 1, 2], [5, 2, 6],  # top face (x-y)
                    [3, 0, 4], [3, 4, 7],  # bottom face (x-y)
                    [4, 0, 5], [1, 5, 0],  # left face (y-z)
                    [2, 3, 7], [2, 7, 6]]  # right face (y-z)

        mesh = Mesh(vertices, triangles)
        self._raysect_mesh = subdivide(mesh)

        
        
