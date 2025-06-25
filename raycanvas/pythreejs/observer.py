import numpy as np

from raycanvas.backend.base.observer import ObserverSource
from raycanvas.pythreejs.visualiser import PythreejsVisualiser

from pythreejs import ArrowHelper, Vector3


class ObserverVisualiser(PythreejsVisualiser):
    """
    Class for graphical representation of the Raysect observers.
    :param raysect_object: Raysect Observer instance
    """

    def __init__(self, source):
        
        self.length = 1e-2
        super().__init__(source)

    def _init_object_kwargs(self):
       self.plot_kwargs = {}
       self.plot_kwargs["headLength"] = 0.2 * self.length 
       self.plot_kwargs["headWidth"] = 0.2 * self.length
       self.plot_kwargs["color"] = "green"
    
    @property
    def length(self):
        return self._length
    
    @length.setter
    def length(self, value):

        if not value > 0:
            raise ValueError("value has to be larger than 0")
        
        self._length = value
        

    def set_source(self, source):

        if not isinstance(source, ObserverSource):
            raise TypeError("source has to be instance of SourceBase.")
        
        self._source = source
    
    def get_data_object(self):
        origin = self._source.origin_root
        origin = (origin.x, origin.y, origin.z)
        direction = self.length * self._source.direction_root.normalise()
        direction= (direction.x, direction.y, direction.z)
        
        return ArrowHelper(dir=direction, origin=origin, length=self.length, **self.plot_kwargs)

    def _add_object(self, plotter):

        mesh = self.get_data_object()

        plotter.scene.add(mesh)

