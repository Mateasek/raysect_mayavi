from raysect_mayavi.primitives.mesh import TriangularMeshSource
from raysect_mayavi.pyvista.visualiser import PyvistaVisualiser

import pyvista as pv
pv.rcParams['use_ipyvtk'] = True


class TriangularMeshVisualiser(PyvistaVisualiser):
    """
    This class serves as base class for Raysect objetcs visualised with triangular meshses.
    """

    def __init__(self, source):

        self._init_plot_kwargs()

        super().__init__(source)

    def set_source(self, source):

        if not isinstance(source, TriangularMeshSource):
            raise TypeError("source has to be instance of TriangularMeshSource.")

        self._source = source

    def get_data_object(self):

        vertices = self._source.vertices
        triangles = self._source.triangles

        return pv.make_tri_mesh(vertices, triangles)

    def _add_object(self, plotter):

        mesh = self.get_data_object()
        plotter.add_mesh(mesh)
