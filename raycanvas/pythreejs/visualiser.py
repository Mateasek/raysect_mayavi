from raysect.core import Point3D

from raycanvas.backend.base.source import SourceBase, VisualiserBase, TriangularMeshSource
from raycanvas.pythreejs.render import get_default_renderer

class PythreejsVisualiser(VisualiserBase):

    """
    This is the base class for the raycanvas representation of raysect objects.
    """
    def __init__(self, source):

        self.set_source(source)
        self._init_figure_kwargs()
        self._init_plot_kwargs()
        self._init_object_kwargs()

    def _init_figure_kwargs(self):
        self.figure_kwargs= {}
        self.figure_kwargs["window_size"] = (512, 256)

    def _init_plot_kwargs(self):
        self.plot_kwargs = {}
    
    def _init_object_kwargs(self):
        self.object_kwargs = {}
        

    @property
    def source(self):
        return self._source
    
    def set_source(self, source):
        raise NotImplementedError("Virtual method _set_source() has not been implemented.")

    def get_data_object(self):
        """
        Constructs the pyvista.DataObject representation of the Raysect primitive
        """
        raise NotImplementedError("Virtual method get_vista_object() has not been implemented.")
    
    def plot(self, plotter=None):
        """
        Plot the Mayavi representation of the Raysect object into the figure. The representation is done always
        in the root node.
        :param figure:Optional, specifies the figure to plot in.
        return mayavi figure
        """
        if plotter is None:
            plotter = get_default_renderer
        #elif not isinstance(plotter, pv.Plotter):
        #    raise ValueError("figure has to be instance of mlab.figure")

        self._add_object(plotter)

        return plotter

    def _add_object(self, plotter):
        raise NotImplementedError("Virtual method has not been implemented.")
