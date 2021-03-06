from raysect_mayavi.primitives.ray import LoggingRaySource

from raysect_mayavi.mayavi.visualiser import MayaviVisualiser

from mayavi import mlab


class LoggingRayVisualiser(MayaviVisualiser):
    """
    Class for graphical representation of the Raysect LoggingRay primitive.
    To display correctly with other scenegraph  components observation has to be usualy
    done in the root node (an instance of World) of the scenegraph.

    :param raysect_object: Raysect Loggingray instance containing path_vertices
    """

    def __init__(self, source):

        self._init_plot_kwargs()

        super().__init__(source)

    def _init_plot_kwargs(self):
        self.plot_kwargs = {}
        self.plot_kwargs["tube_radius"] = None

    def set_source(self, source):

        if not isinstance(source, LoggingRaySource):
            raise TypeError("source has to be instance of LoggingRaySource.")

        self._source = source

    def _mayavi_plot(self, figure):
        vertices = self._source._vertices

        mlab.plot3d(vertices[:, 0], vertices[:, 1], vertices[:, 2],
                    figure=figure, **self.plot_kwargs)
