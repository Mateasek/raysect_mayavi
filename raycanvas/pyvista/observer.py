import numpy as np

from raycanvas.backend.base.observer import ObserverSource
from raycanvas.pyvista.visualiser import PyvistaVisualiser

import pyvista as pv 

class ObserverVisualiser(PyvistaVisualiser):
    """
    Class for graphical representation of the Raysect Box primitive.
    :param raysect_object: Raysect Observer instance
    """

    def __init__(self, source: ObserverSource):
        
        self.length = 0.1
        super().__init__(source)

    def _init_plot_kwargs(self):
       self.plot_kwargs = {}
       self.plot_kwargs["tip_length"] = 2 * self.length 
       self.plot_kwargs["tip_radius"] = 2 * self.length
       self.plot_kwargs["tip_resolution"] = 20
       self.plot_kwargs["shaft_radius"] = 0.2 * self.length
       self.plot_kwargs["shaft_resolution"] = 20
    
    @property
    def length(self):
        return self._length
    
    @length.setter
    def length(self, value):

        if not value > 0:
            raise ValueError("value has to be larger than 0")
        
        self._length = value
        

    def set_source(self, source: ObserverSource):

        if not isinstance(source, ObserverSource):
            raise TypeError("source has to be instance of SourceBase.")
        
        self._source = source
    
    def get_data_object(self) -> pv.Arrow:
        origin_vector = self._source.origin_root
        origin = np.array([origin_vector.x, origin_vector.y, origin_vector.z])
        direction_vector = self.length * self._source.direction_root.normalise()
        direction = np.array([direction_vector.x, direction_vector.y, direction_vector.z])
        
        return pv.Arrow(start=origin, direction=direction, **self.plot_kwargs)

    def _add_object(self, plotter: pv.Plotter):

        mesh = self.get_data_object()
        plotter.add_mesh(mesh, "r")

