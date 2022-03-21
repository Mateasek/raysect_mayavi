
from raysect_mayavi.primitives.source import SourceBase, TriangularMeshSource
from raysect_mayavi.primitives.ray import LoggingRaySource
from raysect_mayavi.primitives.observer import ObserverSource

from raysect_mayavi.pyvista.mesh import TriangularMeshVisualiser
from raysect_mayavi.pyvista.ray import LoggingRayVisualiser
from raysect_mayavi.pyvista.observer import ObserverVisualiser


def parse_sources(sources):
    
    visualisers = {}
    
    for key, source in sources.items():
        visualisers[key] = assign_visualiser(source)

    
    return visualisers

def assign_visualiser(source):
    
    if not isinstance(source, SourceBase):
        raise TypeError("sources items has to be of type SourceBase.")
    if isinstance(source, TriangularMeshSource):
        return TriangularMeshVisualiser(source)
    elif isinstance(source, LoggingRaySource):
        return LoggingRayVisualiser(source)
    elif isinstance(source, ObserverSource):
        return ObserverVisualiser(source)
    else:
        raise TypeError("Source type not recognised")
    