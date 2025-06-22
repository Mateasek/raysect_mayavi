import numpy as np

from raysect_mayavi.primitives.ray import LoggingRaySource
from raysect_mayavi.pythreejs.visualiser import PythreejsVisualiser

from pythreejs import (LineSegmentsGeometry, LineSegments2, LineMaterial, BufferGeometry,
                       BufferAttribute)


class LoggingRayVisualiser(PythreejsVisualiser):
    """
    Class for graphical representation of the Raysect LoggingRay primitive. To display correctly with other scenegraph 
    components observation has to be usualy done in the root node (an instance of World) of the scenegraph.
    :param raysect_object: Raysect Loggingray instance containing path_vertices
    """

    def __init__(self, source):
        
        self._init_plot_kwargs()
        super().__init__(source)

    def _init_object_kwargs(self):
        self.object_kwargs = {}
        self.object_kwargs["color"] = "red"
        self.object_kwargs["linewidth"] = 1
        self.continuum_length = 0.1

    @property
    def continuum_length(self):
        return self._contunuum_length

    @continuum_length.setter
    def continuum_length(self, value):
        if not value > 0:
            raise ValueError("Value has to be larger than 0.")
        
        self._contunuum_length = value

    def set_source(self, source):

        if not isinstance(source, LoggingRaySource):
            raise TypeError("source has to be instance of SourceBase.")
        
        self._source = source

    def get_data_object(self):
        vertices = self._source._vertices.astype(np.float32)

        if len(vertices) == 0:
            origin = self._source._raysect_object.origin
            direction = self._source._raysect_object.direction
            ray_end = origin + direction.normalise() * self._contunuum_length

            vertices = np.array([[*origin], [*ray_end]], dtype=np.float32)
        
        line_segments = np.zeros((vertices.shape[0]-1, 2, 3), dtype=np.float32)

        for i in range(line_segments.shape[0]):
            line_segments[i, 0, :] = vertices[i, :]
            line_segments[i, 1, :] = vertices[i+1, :]

        geometry = LineSegmentsGeometry(positions=line_segments)
        material = LineMaterial(**self.object_kwargs)

        ray_path = LineSegments2(geometry, material)

        return ray_path

    def _add_object(self, plotter):
        
        mesh = self.get_data_object()
        if mesh is not None:
            return plotter.scene.add(mesh)