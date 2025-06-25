from raycanvas.backend.base.source import SourceBase, TriangularMeshSource
from raycanvas.backend.base.ray import LoggingRaySource
from raycanvas.backend.base.observer import ObserverSource
from raycanvas.pyvista.visualiser import PyvistaVisualiser
from raycanvas.pyvista.mesh import TriangularMeshVisualiser
from raycanvas.pyvista.ray import LoggingRayVisualiser
from raycanvas.pyvista.observer import ObserverVisualiser


def parse_sources(sources: dict) -> dict:
    
    visualisers = {}
    
    for key, source in sources.items():
        visualisers[key] = assign_visualiser(source)

    
    return visualisers

def assign_visualiser(source: SourceBase) -> PyvistaVisualiser:
    
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
    