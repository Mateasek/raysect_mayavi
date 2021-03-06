from raysect_mayavi.primitives.mesh import TriangularMeshSource
from raysect_mayavi.mayavi.visualiser import MayaviVisualiser

from mayavi import mlab


class TriangularMeshVisualiser(MayaviVisualiser):
    """
    This class serves as base class for Raysect objetcs visualised with triangular meshses.
    """

    def __init__(self, source):

        self._init_plot_kwargs()

        super().__init__(source)
        self.plot_method = mlab.triangular_mesh

    def set_source(self, source):

        if not isinstance(source, TriangularMeshSource):
            raise TypeError("source has to be instance of TriangularMeshSource.")

        self._source = source

    def _init_plot_kwargs(self):
        self.plot_kwargs = {}
        self.plot_kwargs["color"] = (163 / 255.0, 163 / 255.0, 163 / 255.0)
        self.plot_kwargs["transparent"] = True
        self.plot_kwargs["opacity"] = 1

    def _mayavi_plot(self, figure):
        vertices = self._source.vertices

        mlab.triangular_mesh(vertices[:, 0], vertices[:, 1], vertices[:, 2], self._source.triangles,
                             figure=figure, **self.plot_kwargs)
