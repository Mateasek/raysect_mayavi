from raysect.core import Point3D
from raysect.primitive import Mesh

from mayavi import mlab
from mayavi.core.module import Module as MayaviModule
from mayavi.core.scene import Scene

class SourceBase:
    """
    This is the base class for representation of raysect objects in visualisation toolkits.
    """
    def __init__(self, raysect_object):

        self._raysect_object = raysect_object
        self._graphic_source_from_raysect_object()

    def _graphic_source_from_raysect_object(self):
        """
        Constructs source of information for visualisation from raysect primitives.
        :return:
        """
        raise NotImplementedError("Virtual method _graphic_source_from_raysect_object() has not been implemented.")

class VisualiserBase:
    
    def __init__(self, primitive_source):
        
        if not isinstance(primitive_source, SourceBase):
            raise TypeError("primitive_sourse has to be instance of SourceBase")

        self._primitive_source = primitive_source
        
    @property
    def primitive_source(self):
        return self._primitive_source
    
    def plot(self):
        raise NotImplementedError("Virtual method plot() has not been implemented.")
        

class TriangularMeshSource(SourceBase):
    """
    This class serves as base class for Raysect objetcs visualised with triangular meshses.
    """
    def __init__(self, raysect_object):

        self._raysect_mesh = None
        super().__init__(raysect_object)

    def _vertices_root(self): 
        if self._raysect_object.parent:
            to_world = self._raysect_object.to_root()
        else:
            to_world = self._raysect_object.transform
        return self._vertices_transform(to_world)

    def vertices_node(self, node):
        transform = self._raysect_object.to(node) 
        return self._vertices_transform(transform)

    def _vertices_transform(self, transform):
        vertices = self._raysect_mesh.data.vertices.copy()

        for i in range(vertices.shape[0]):
            point = Point3D(vertices[i, 0], vertices[i, 1], vertices[i, 2]).transform(transform)
            vertices[i, 0] = point.x
            vertices[i, 1] = point.y
            vertices[i, 2] = point.z

        return vertices

    @property
    def raysect_object(self):
        return self._raysect_object

    @property
    def vertices_local(self):
        return self._raysect_mesh.data.vertices()

    @property
    def vertices(self):
        return self._vertices_root()

    @property
    def vertices_x(self):
        return self._vertices_root()[:, 0]

    @property
    def vertices_y(self):
        return self._vertices_root()[:, 1]

    @property
    def vertices_z(self):
        return self._vertices_root()[:, 2]

    @property
    def triangles(self):
        return self._raysect_mesh.data.triangles[:,0:3]
