"""Module for visualizing Raysect scenegraph using PyVista.

This module provides functions to visualize Raysect scenegraphs and objects
using PyVista's visualization capabilities.
"""

import pyvista as pv

from raysect.core import Node

from raycanvas.backend.builtin.primitives.parse_nodes import parse_nodes, parse_elements
from raycanvas.pyvista.parse import parse_sources


def get_scenegraph_visualisers(node: Node) -> dict:
    """Extract visualisers from a Raysect scenegraph.

    Args:
        node: A Raysect Node instance.

    Returns:
        dict: A dictionary of visualisers for the nodes in the node.
    """
    sources = parse_nodes(node)
    return parse_sources(sources)


def get_visualisers(object_list: list) -> dict:
    """Extract visualisers from a list of Raysect objects.

    Args:
        object_list: A list of Raysect objects to visualize.

    Returns:
        dict: A dictionary of visualisers for the objects.
    """
    sources = parse_elements(object_list)
    return parse_sources(sources)


def visualise(
    object_list: list,
    plotter: pv.Plotter = None,
    show_axes: bool = False,
    axes_length: float = 1,
    jupyter_backend: str = "trame",
) -> pv.Plotter:
    """Visualize a list of Raysect objects.

    Args:
        object_list: A list of Raysect objects to visualize.
        plotter: Optional PyVista plotter instance. If None, a new one will be created.
        show_axes: Whether to show coordinate axes (default: False).
        axes_length: Length of the coordinate axes (default: 1).
        jupyter_backend: Backend to use for Jupyter visualization (default: "trame").

    Returns:
        pv.Plotter: The PyVista plotter instance containing the visualization.
    """
    pv.set_jupyter_backend(jupyter_backend)
    visualisers = get_scenegraph_visualisers(object_list)
    return visualise_visualisers(visualisers, plotter, show_axes, axes_length)


def visualise_scenegraph(
    node: Node,
    plotter: pv.Plotter = None,
    show_axes: bool = False,
    axes_length: float = 1,
    jupyter_backend: str = "trame",
) -> pv.Plotter:
    """Visualize a Raysect scenegraph.

    Args:
        node: A Raysect Node instance.
        plotter: Optional PyVista plotter instance. If None, a new one will be created.
        show_axes: Whether to show coordinate axes (default: False).
        axes_length: Length of the coordinate axes (default: 1).
        jupyter_backend: Backend to use for Jupyter visualization (default: "trame").

    Returns:
        pv.Plotter: The PyVista plotter instance containing the visualization.
    """
    pv.set_jupyter_backend(jupyter_backend)
    visualisers = get_scenegraph_visualisers(node)
    return visualise_visualisers(visualisers, plotter, show_axes, axes_length)


def visualise_visualisers(
    visualisers: dict,
    plotter: pv.Plotter = None,
    show_axes: bool = False,
    axes_length: float = 1,
) -> pv.Plotter | None:
    """Visualize a dictionary of visualisers.

    Args:
        visualisers: Dictionary of visualisers to plot.
        plotter: Optional PyVista plotter instance. If None, a new one will be created.
        show_axes: Whether to show coordinate axes (default: False).
        axes_length: Length of the coordinate axes (default: 1).

    Returns:
        pv.Plotter: The PyVista plotter instance containing the visualization.
    """
    if plotter is None:
        plotter = pv.Plotter(window_size=(1024, 768), notebook=True, off_screen=False)

    if show_axes:
        # Add coordinate axes
        plotter.add_mesh(pv.Arrow((0, 0, 0), (axes_length, 0, 0)), color="r")  # X-axis
        plotter.add_mesh(pv.Arrow((0, 0, 0), (0, axes_length, 0)), color="g")  # Y-axis
        plotter.add_mesh(pv.Arrow((0, 0, 0), (0, 0, axes_length)), color="b")  # Z-axis

    # Plot each visualiser
    for visualiser in visualisers.values():
        visualiser.plot(plotter)

    # Display the plotter in Jupyter notebook
    plotter.show()

    return plotter
