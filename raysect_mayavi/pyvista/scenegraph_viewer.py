import pyvista as pv

from raysect_mayavi.primitives.parse_nodes import parse_nodes
from raysect_mayavi.pyvista.parse import parse_sources


def visualise_scenegraph(world, plotter=None, show_axes=False, axes_length=1,
                         jupyter_backend="pythreejs"):
    """
    Automatic visualisation of the raysect scenegraph.
    :param world: Instance of the Raysect World 
    """

    pv.set_jupyter_backend(jupyter_backend)

    if plotter is None:
        plotter = pv.Plotter(window_size=(1024, 768))

    sources = parse_nodes(world)
    visualisers = parse_sources(sources)

    for _, visualiser in visualisers.items():

        visualiser.plot(plotter)


    return plotter