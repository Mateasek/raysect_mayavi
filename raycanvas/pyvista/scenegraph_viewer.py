from typing import Union

import pyvista as pv

from raysect.core import Node, Primitive, Observer
from raysect.optical.loggingray import LoggingRay

from raycanvas.backend.builtin.primitives.parse_nodes import parse_elements
from raycanvas.pyvista.parse import parse_sources
from raycanvas.pyvista.visualiser import PyvistaVisualiser

_uni_type_raysect_elements = Union[Node, list[Node], Primitive, Observer, LoggingRay]


def get_visualisers(elements: _uni_type_raysect_elements) -> dict:
    """Extract visualisers from a list of Raysect objects.

    Args:
        object_list: A list of Raysect objects to visualize.

    Returns:
        dict: A dictionary of visualisers for the objects.
    """
    sources = parse_elements(elements)
    return parse_sources(sources)


def visualise(
    elements: _uni_type_raysect_elements,
    visualisers: dict[Node, PyvistaVisualiser] | None = None,
    plotter: pv.Plotter | None = None,
    show_axes: bool = False,
    axes_length: float = 1,
    jupyter_backend: str = "trame",
) -> tuple[pv.Plotter, dict[Node, PyvistaVisualiser]]:
    """Visualize elements of a Raysect scenegraph.

    Args:
        elements: Raysect scenegraph elements to visualize.
        visualisers: Dictionary of raycanvas visualisers to plot. If not None, new visualisers will be added to the existing dictionary.
        plotter: Optional PyVista plotter instance. If None, a new one will be created.
        show_axes: Whether to show coordinate axes (default: False).
        axes_length: Length of the coordinate axes (default: 1).

    Returns:
        pv.Plotter: The PyVista plotter instance containing the visualization.
    """

    new_visualisers = get_visualisers(elements)

    if visualisers is None:
        visualisers = new_visualisers
    else:
        visualisers = {**visualisers, **new_visualisers}

    if plotter is None:
        plotter = pv.Plotter(window_size=(1024, 768), notebook=True, off_screen=False)
        pv.set_jupyter_backend(jupyter_backend)

        if show_axes:
            # Add coordinate axes
            plotter.add_mesh(
                pv.Arrow((0, 0, 0), (axes_length, 0, 0)), color="r"
            )  # X-axis
            plotter.add_mesh(
                pv.Arrow((0, 0, 0), (0, axes_length, 0)), color="g"
            )  # Y-axis
            plotter.add_mesh(
                pv.Arrow((0, 0, 0), (0, 0, axes_length)), color="b"
            )  # Z-axis

    # Plot each visualiser
    for new_visualiser in new_visualisers.values():
        new_visualiser.plot(plotter)

    return plotter, visualisers
