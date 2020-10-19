
import numpy as np

from raysect.core import Point3D, Vector3D
from raysect.core.scenegraph.observer import Observer
from raysect_mayavi.primitives.source import MayaviSource

from mayavi import mlab


class ObserverSource(MayaviSource):
    """
    Representation of a Raysect Observer objects.
    :param raysect_object: Raysect Observer
    """

    def __init__(self, raysect_object):
        
        if not isinstance(raysect_object, Observer):
            raise TypeError("The raysect_object has to be instance of Raysect Observer primitive, wrong type '{}' given.".format(type(raysect_object)))

        self._init_plot_kwargs()
        
        super().__init__(raysect_object)
        
    def _init_plot_kwargs(self):
        self.plot_kwargs = {}
        self.plot_kwargs["scale_factor"] = 0.1

    def _mayavi_source_from_raysect_object(self):
        pass

    def _mayavi_plot(self, figure):
        origin_raysect = Point3D(0, 0, 0).transform(self._raysect_object.to_root())
        direction_raysect = Vector3D(0, 0, 1).transform(self._raysect_object.to_root())

        x = np.array([origin_raysect.x], ndmin=1)
        y = np.array([origin_raysect.y], ndmin=1)
        z = np.array([origin_raysect.z], ndmin=1)

        u = np.array([direction_raysect.x], ndmin=1)
        v = np.array([direction_raysect.y], ndmin=1)
        w = np.array([direction_raysect.z], ndmin=1)


        mlab.quiver3d(x, y, z, u, v, w,
                    figure=figure, **self.plot_kwargs)
