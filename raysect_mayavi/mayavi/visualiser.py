from raysect.core import Point3D

from raysect_mayavi.primitives.source import SourceBase, VisualiserBase, TriangularMeshSource

from mayavi import mlab
from mayavi.core.scene import Scene

class MayaviVisualiser(VisualiserBase):

    """
    This is the base class for the raysect_mayavi representation of raysect objects.
    """
    def __init__(self, source):

        self.set_source(source)
        self._init_figure_kwargs()

    def _init_figure_kwargs(self):
        self.figure_kwargs= {}
        self.figure_kwargs["size"] = (1024, 768)
        self.figure_kwargs["bgcolor"] = (1, 1, 1)
        self.figure_kwargs["fgcolor"] = (0.5, 0.5, 0.5)

    @property
    def source(self):
        return self._source
    
    def set_source(self, source):
        raise NotImplementedError("Virtual method _set_source() has not been implemented.")
    
    def plot(self, figure=None):
        """
        Plot the Mayavi representation of the Raysect object into the figure. The representation is done always
        in the root node.
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
