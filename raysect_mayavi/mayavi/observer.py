import numpy as np

from raysect_mayavi.primitives.observer import ObserverSource
from raysect_mayavi.mayavi.visualiser import MayaviVisualiser

from mayavi import mlab


class ObserverVisualiser(MayaviVisualiser):
    """
    Class for graphical representation of the Raysect Box primitive.

    :param raysect_object: Raysect Observer instance
    """

    def __init__(self, source):

        self._init_plot_kwargs()

        super().__init__(source)

    def _init_plot_kwargs(self):
        self.plot_kwargs = {}
        self.plot_kwargs["scale_factor"] = 0.1

    def set_source(self, source):

        if not isinstance(source, ObserverSource):
            raise TypeError("source has to be instance of ObserverSource.")

        self._source = source

    def _mayavi_plot(self, figure):
        origin_raysect = self._source.origin_root
        direction_raysect = self._source.direction_root

        x = np.array([origin_raysect.x], ndmin=1)
        y = np.array([origin_raysect.y], ndmin=1)
        z = np.array([origin_raysect.z], ndmin=1)

        u = np.array([direction_raysect.x], ndmin=1)
        v = np.array([direction_raysect.y], ndmin=1)
        w = np.array([direction_raysect.z], ndmin=1)

        mlab.quiver3d(x, y, z, u, v, w,
                      figure=figure, **self.plot_kwargs)
