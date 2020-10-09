
from mayavi import mlab

from raysect.core import Observer, World, Point3D, Vector3D, Node, Primitive
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
        
        if isinstance(child, Primitive):
            mesh_dict[child] = to_mesh(child)

    return mesh_dict


def visualise_scenegraph(world, figure=None, show_axes=False, axes_length=1):

    if not isinstance(world, World):
        raise TypeError("The visualisation function takes a Raysect World object as its argument.")
    
    if figure is None:
        figure = mlab.figure(size=(1024, 768), bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5))

    meshes = parse_nodes(world)

    for raysect_primitive, mesh_description in meshes.items():

        mesh_description.mayavi_plot(figure)

    if show_axes:
        mlab.plot3d([0, axes_length], [0, 0], [0, 0], tube_radius=axes_length/100, color=(1, 0, 0))
        mlab.plot3d([0, 0], [0, axes_length], [0, 0], tube_radius=axes_length/100, color=(0, 1, 0))
        mlab.plot3d([0, 0], [0, 0], [0, axes_length], tube_radius=axes_length/100, color=(0, 0, 1))

    return figure