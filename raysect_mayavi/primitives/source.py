from raysect.core import Point3D
from raysect.primitive import Mesh

from mayavi import mlab
from mayavi.core.module import Module as MayaviModule
from mayavi.core.scene import Scene

class MayaviSource:
    """
    This is the base class for the raysect_mayavi representation of raysect objects.
    """
    def __init__(self, raysect_object):

        self._init_figure_kwargs()
        self._raysect_object = raysect_object
        self._mayavi_source_from_raysect_object()

    def _init_figure_kwargs(self):
        self.figure_kwargs= {}
        self.figure_kwargs["size"] = (1024, 768)
        self.figure_kwargs["bgcolor"] = (1, 1, 1)
        self.figure_kwargs["fgcolor"] = (0.5, 0.5, 0.5)
        

    def _mayavi_source_from_raysect_object(self):
        """
        Constructs raysect Mesh from sources.
        :return:
        """
        raise NotImplementedError("Virtual method _mayavi_source_from_raysect_object() has not been implemented.")

    def mayavi_plot(self, figure=None):
        """
        Plot the Mayavi representation of the Raysect object into the figure. The representation is done always in the root node.
        :param figure:Optional, specifies the figure to plot in.
        return mayavi figure
        """
        if figure is None:
            figure = mlab.figure(**self.figure_kwargs)
        elif not isinstance(figure, Scene):
            raise ValueError("figure has to be instance of mlab.figure")

        self._mayavi_plot(figure)

        return figure
    
    def _mayavi_plot(self, figure):
        raise NotImplementedError("Virtual method _mayavi_plot() has not been implemented.")

class TriangularMeshSource(MayaviSource):

    def __init__(self, raysect_object):

        self._raysect_mesh = None
        self._init_plot_kwargs()

        super().__init__(raysect_object)
        
        self.plot_method = mlab.triangular_mesh    
         
    def _init_plot_kwargs(self):
        self.plot_kwargs = {}
        self.plot_kwargs["color"] = (163 / 255.0, 163 / 255.0, 163 / 255.0)
        self.plot_kwargs["transparent"] = True
        self.plot_kwargs["opacity"] = 1

    def _vertices_root(self): 
        if self._raysect_object.parent:
            to_world = self._raysect_object.to_root()
        else:
            to_world = self._raysect_object.transform
        return self._vertices_transform(to_world)
    
    def vertices_node(self, node):
        transform = self._raysect_object.to(node) 
        return self._vertices_transform(to_world)
     
    def _vertices_transform(self, transform):
        vertices = self._raysect_mesh.data.vertices.copy()


        for i in range(vertices.shape[0]):
            p = Point3D(vertices[i, 0], vertices[i, 1], vertices[i, 2]).transform(transform)
            vertices[i, 0] = p.x
            vertices[i, 1] = p.y
            vertices[i, 2] = p.z

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

    def _mayavi_plot(self, figure):
        vertices = self._vertices_root()

        mlab.triangular_mesh(vertices[:,0], vertices[:,1], vertices[:,2], self.triangles,
                             figure=figure, **self.plot_kwargs)