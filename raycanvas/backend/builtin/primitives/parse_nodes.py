from typing import Union
from raysect.core import Observer, World, Node, Primitive
from raysect.optical.loggingray import LoggingRay

from raycanvas.backend.builtin.primitives.mesh import to_mesh, MeshSource

# Use consistent type annotation syntax
RaysectElement = Union[Node, Primitive, Observer, LoggingRay]
RaysectElements = Union[RaysectElement, list[RaysectElement]]


def parse_elements(
    element: RaysectElements, mesh_dict: dict[Node, MeshSource] | None = None
) -> dict[Node, MeshSource]:
    """Parse Raysect scenegraph elements and create mesh representations.

    Searches through Raysect scenegraph elements and creates a dictionary with
    primitive representations for visualization.

    Args:
        element: A Raysect element or list of elements to be parsed
        mesh_dict: Optional dictionary to add primitive representations to

    Returns:
        Dictionary with keys being raysect primitives and values being their mesh representations

    Raises:
        TypeError: If element is not a valid Raysect type
    """
    if mesh_dict is None:
        mesh_dict = {}

    if isinstance(element, list):
        # Process each element in the list
        for obj in element:
            mesh_dict = parse_raysect_object(obj, mesh_dict)
    else:
        mesh_dict = parse_raysect_object(element, mesh_dict)

    return mesh_dict


def parse_raysect_object(
    element: RaysectElement, mesh_dict: dict[Node, MeshSource] | None = None
) -> dict[Node, MeshSource]:
    """Parse a single Raysect object and add its representation.

    Parses through an underlying Raysect scenegraph branch if element is of type Node.

    Args:
        element: A single Raysect element to be parsed
        mesh_dict: Optional dictionary to add primitive representations to

    Returns:
        Dictionary with keys being raysect primitives and values being their mesh representations

    Raises:
        TypeError: If element is not a valid Raysect type
    """
    if not isinstance(element, (Node, World, Primitive, Observer, LoggingRay)):
        raise TypeError(
            "element must be an instance of raysect.core.Node, raysect.core.World, "
            "raysect.core.Primitive, or raysect.core.Observer"
        )

    if mesh_dict is None:
        mesh_dict = {}

    # Handle Node types by recursively processing children
    if isinstance(element, (Node, World)):
        for child in element.children:
            mesh_dict = parse_raysect_object(child, mesh_dict)

    # Handle primitive types by creating mesh representations
    if isinstance(element, (Primitive, Observer, LoggingRay)):
        mesh_dict[element] = to_mesh(element)

    return mesh_dict
