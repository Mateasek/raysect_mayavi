
from mayavi import mlab

from raysect.core import Observer, World, Point3D, Vector3D, Node
from raysect_mayavi.primitives import to_mesh


def _parse_nodes(node, mesh_list):

    for child in node.children:

        if type(child) == Node:
            _parse_nodes(child, mesh_list)
        elif isinstance(child, Observer):
            return
        else:
            mesh_list.append((to_mesh(child), child.meta))


def visualise_scenegraph(world, show_axes=False, axes_length=1):

    if not isinstance(world, World):
        raise TypeError("The visualisation function takes a Raysect World object as its argument.")

    fig = mlab.figure(size=(1024, 768), bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5))

    meshes = []
    _parse_nodes(world, meshes)

    for mesh_description in meshes:

        vertices, triangles = mesh_description[0]
        dx = vertices[:, 0]
        dy = vertices[:, 1]
        dz = vertices[:, 2]

        meta = mesh_description[1]
        try:
            color = meta['viz-color']
        except KeyError:
            color = (163 / 255.0, 163 / 255.0, 163 / 255.0)
        try:
            opacity = meta['viz-opacity']
            transparent = True
        except KeyError:
            opacity = 1
            transparent = False

        mlab.triangular_mesh(dx, dy, dz, triangles, color=color,
                             figure=fig, transparent=transparent, opacity=opacity)

    if show_axes:
        mlab.plot3d([0, axes_length], [0, 0], [0, 0], tube_radius=axes_length/100, color=(1, 0, 0))
        mlab.plot3d([0, 0], [0, axes_length], [0, 0], tube_radius=axes_length/100, color=(0, 1, 0))
        mlab.plot3d([0, 0], [0, 0], [0, axes_length], tube_radius=axes_length/100, color=(0, 0, 1))
