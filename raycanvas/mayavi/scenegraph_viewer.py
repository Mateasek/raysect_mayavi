from mayavi import mlab

from raycanvas.backend.builtin.primitives.parse_nodes import parse_nodes
from raycanvas.mayavi.parse import parse_sources


def visualise_scenegraph(world, figure=None, show_axes=False, axes_length=1):
    """
    Automatic visualisation of the raysect scenegraph.
    :param world: Instance of the Raysect World 
    """

    if figure is None:
        figure = mlab.figure(size=(1024, 768), bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5))

    sources = parse_nodes(world)
    visualisers = parse_sources(sources)

    for _, visualiser in visualisers.items():

        visualiser.plot(figure)

    if show_axes:
        mlab.plot3d([0, axes_length], [0, 0], [0, 0], tube_radius=axes_length/100, color=(1, 0, 0))
        mlab.plot3d([0, 0], [0, axes_length], [0, 0], tube_radius=axes_length/100, color=(0, 1, 0))
        mlab.plot3d([0, 0], [0, 0], [0, axes_length], tube_radius=axes_length/100, color=(0, 0, 1))

    return figure