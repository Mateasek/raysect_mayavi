import numpy as np

from raycanvas.backend.builtin.primitives.mesh import TriangularMeshSource
from raycanvas.pyvista.visualiser import PyvistaVisualiser

import pyvista as pv


class TriangularMeshVisualiser(PyvistaVisualiser):
    """
    This class serves as base class for Raysect objetcs visualised with triangular meshses.
    """
    def __init__(self, source: TriangularMeshSource):

        
        self._init_plot_kwargs()

        
        super().__init__(source)

    def set_source(self, source: TriangularMeshSource):

        if not isinstance(source, TriangularMeshSource):
            raise TypeError("source has to be instance of SourceBase.")
        
        self._source = source
    
    def get_data_object(self) -> pv.PolyData:

        vertices = self._source.vertices
        triangles = self._source.triangles

        return pv.make_tri_mesh(vertices, triangles)

        
    def _add_object(self, plotter: pv.Plotter):
        
        mesh = self.get_data_object()
        plotter.add_mesh(mesh, color=None)