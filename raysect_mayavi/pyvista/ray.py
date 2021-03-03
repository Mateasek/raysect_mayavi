import numpy as np

from raysect_mayavi.primitives.ray import LoggingRaySource
from raysect_mayavi.pyvista.visualiser import PyvistaVisualiser

import pyvista as pv 
pv.rcParams['use_ipyvtk'] = True

class LoggingRayVisualiser(PyvistaVisualiser):
    """
    Class for graphical representation of the Raysect LoggingRay primitive. To display correctly with other scenegraph 
    components observation has to be usualy done in the root node (an instance of World) of the scenegraph.
    :param raysect_object: Raysect Loggingray instance containing path_vertices
    """

    def __init__(self, source):
        
        self._init_plot_kwargs()
        super().__init__(source)

    def set_source(self, source):

        if not isinstance(source, LoggingRaySource):
            raise TypeError("source has to be instance of SourceBase.")
        
        self._source = source

    def get_data_object(self):
        vertices = self._source._vertices
        if len(vertices) > 0:
            return pv.lines_from_points(vertices, close=False)
        else:
            return None

    def _add_object(self, plotter):

        mesh = self.get_data_object()
        if mesh is not None:
            plotter.add_mesh(mesh, color="g")