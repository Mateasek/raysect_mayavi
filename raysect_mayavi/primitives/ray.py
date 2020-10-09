import numpy as np

from raysect.optical.loggingray import LoggingRay
from raysect_mayavi.primitives.source import MayaviSource

from mayavi import mlab


class LoggingRaySource(MayaviSource):
    """
    Representation of the Raysect LoggingRay.
    :param raysect_object: Raysect LoggingRay
    """

    def __init__(self, raysect_object):
        
        if not isinstance(raysect_object, LoggingRay):
            raise TypeError("The raysect_object has to be instance of Raysect Box primitive, wrong type '{}' given.".format(type(raysect_object)))

        self._init_plot_kwargs()
        
        super().__init__(raysect_object)
        
    def _init_plot_kwargs(self):
        self.plot_kwargs = {}
        self.plot_kwargs["tube_radius"] = None

    def _mayavi_source_from_raysect_object(self):

        self._vertices = np.zeros((len(self._raysect_object.path_vertices), 3))
        for i, vertex_point in enumerate(self._raysect_object.path_vertices):
            self._vertices[i, 0] = vertex_point.x
            self._vertices[i, 1] = vertex_point.y
            self._vertices[i, 2] = vertex_point.z

    def _mayavi_plot(self, figure):
        mlab.plot3d(self._vertices[:,0], self._vertices[:,1], self._vertices[:,2],
                    figure=figure, **self.plot_kwargs)

        