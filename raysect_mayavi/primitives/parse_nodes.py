
from raysect.core import Observer, World, Node, Primitive
from raysect_mayavi.primitives import to_mesh


def parse_nodes(node, mesh_dict=None):
    """
    The method searches a Raysect scenegraph and creates a dictionary with primitive representations.
    :param node: A raysect.core.World or raysect.core.Node instance to be searched
    :param mesh_dict: A python dictionary to add primitive representations to
    
    :return A python dictionary with keys being raysect primitives and items being their representations.
    """

    if not isinstance(node, Node) and not isinstance(node, World):
        raise TypeError("node has to be an instance of raysect.core.Node or raysect.core.World")

    if mesh_dict is None:
        mesh_dict = {}

    for child in node.children:
        mesh_dict = parse_nodes(child, mesh_dict)

        if isinstance(child, Primitive) or isinstance(child, Observer):
            mesh_dict[child] = to_mesh(child)

    return mesh_dict