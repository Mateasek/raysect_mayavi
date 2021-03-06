from raysect_mayavi.primitives.parse_nodes import parse_nodes
from raysect_mayavi.pyvista.parse import parse_sources

import pyvista as pv
pv.rcParams['use_ipyvtk'] = True


def visualise_scenegraph(world, plotter=None, show_axes=False, axes_length=1):
    """
    Automatic visualisation of the raysect scenegraph.

    :param world: Instance of the Raysect World
    :param plotter: Optional, specifies the pyvista.plotter to plot in.

    :return: pyvista.plotter
    """

    plotter = plotter or pv.Plotter(window_size=(1024, 768))
    if not isinstance(plotter, pv.Plotter):
        raise ValueError("plotter has to be instance of pyvista.Plotter.")

    sources = parse_nodes(world)
    visualisers = parse_sources(sources)

    for _, visualiser in visualisers.items():

        visualiser.plot(plotter)

    return plotter
