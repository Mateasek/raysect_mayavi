import numpy as np

from raycanvas.backend.base.ray import LoggingRaySource
from raycanvas.pyvista.visualiser import PyvistaVisualiser

import pyvista as pv 

class LoggingRayVisualiser(PyvistaVisualiser):
    """
    Class for graphical representation of the Raysect LoggingRay primitive. To display correctly with other scenegraph 
    components observation has to be usualy done in the root node (an instance of World) of the scenegraph.
    :param raysect_object: Raysect Loggingray instance containing path_vertices
    """

    def __init__(self, source: LoggingRaySource):
        
        self._init_plot_kwargs()
        super().__init__(source)

    def set_source(self, source: LoggingRaySource):

        if not isinstance(source, LoggingRaySource):
            raise TypeError("source has to be instance of SourceBase.")
        
        self._source = source

    def get_data_object(self) -> pv.Line:
        vertices = self._source._vertices

        if len(vertices) < 1:
            return None
        
        # this is a hack, for some reason pythreejs doesnt display lines...
        #lines = np.zeros((vertices.shape[0] - 1, 3), dtype=int)
        #for i in range(vertices.shape[0] - 1):
        #    lines[i, 0] = i
        #    lines[i, 1] = i + 1
        #    lines[i, 2] = i

        #mesh = pv.make_tri_mesh(vertices, lines)
        mesh = pv.MultipleLines(points=vertices)

        return mesh

    def _add_object(self, plotter):

        path = self.get_data_object()
        if path is not None:
            plotter.add_mesh(path, **self.plot_kwargs)