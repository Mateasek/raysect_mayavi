import numpy as np

from raysect.core import Point3D, Vector3D
from raysect.core.scenegraph.observer import Observer
from raysect_mayavi.primitives.source import SourceBase

class ObserverSource(SourceBase):
    """
    Class for graphical representation of the Raysect Box primitive.
    :param raysect_object: Raysect Observer instance
    """

    def __init__(self, raysect_object):
        
        if not isinstance(raysect_object, Observer):
            raise TypeError("The raysect_object has to be instance of Raysect Observer primitive, wrong type '{}' given.".format(type(raysect_object)))

        super().__init__(raysect_object)

    def _graphic_source_from_raysect_object(self):
        pass

    @property
    def origin_root(self):
        """Point3D(0, 0, 1) transformed to the root node."""
        return Point3D(0, 0, 0).transform(self._raysect_object.to_root())

    @property
    def direction_root(self):
        """Vector(0, 0, 1) transformed to the root node."""
        return Vector3D(0, 0, 1).transform(self._raysect_object.to_root())
