import numpy as np

from raysect.optical.loggingray import LoggingRay
from raysect_mayavi.primitives.source import SourceBase


class LoggingRaySource(SourceBase):
    """
    Class for graphical representation of the Raysect LoggingRay primitive. To display correctly with other scenegraph 
    components observation has to be usualy done in the root node (an instance of World) of the scenegraph.
    :param raysect_object: Raysect Loggingray instance containing path_vertices
    """

    def __init__(self, raysect_object):

        if not isinstance(raysect_object, LoggingRay):
            raise TypeError("The raysect_object has to be instance of LoggingRay"
                            " primitive, wrong type '{}' given.".format(type(raysect_object)))

        super().__init__(raysect_object)

    def _graphic_source_from_raysect_object(self):

        self._vertices = np.zeros((len(self._raysect_object.path_vertices), 3))
        for i, vertex_point in enumerate(self._raysect_object.path_vertices):
            self._vertices[i, 0] = vertex_point.x
            self._vertices[i, 1] = vertex_point.y
            self._vertices[i, 2] = vertex_point.z