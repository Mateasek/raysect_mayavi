import pyvista as pv
pv.rcParams['use_ipyvtk'] = True

from raysect_mayavi.primitives.parse_nodes import parse_nodes
from raysect_mayavi.mayavi.parse import parse_sources


def visualise_scenegraph(world, plotter=None, show_axes=False, axes_length=1):
    """
    Automatic visualisation of the raysect scenegraph.
    :param world: Instance of the Raysect World 
    """

    if plotter is None:
        figure = pv.Plotter(window_size=(1024, 768))

    sources = parse_nodes(world)
    visualisers = parse_sources(sources)

    for _, visualiser in visualisers.items():

        visualiser.plot(figure)


    return plotter