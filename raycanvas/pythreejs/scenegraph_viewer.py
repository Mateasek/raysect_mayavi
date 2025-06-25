from raycanvas.backend.builtin.primitives.parse_nodes import parse_nodes
from raycanvas.pythreejs.parse import parse_sources
from raycanvas.pythreejs.render import get_default_camera, get_default_renderer, get_default_scene
from pythreejs import Renderer, OrbitControls, PerspectiveCamera, Scene, AxesHelper



def visualise_scenegraph(world, renderer=None, camera=None, scene=None, axes=True):
    """
    Automatic visualisation of the raysect scenegraph.
    :param world: Instance of the Raysect World 
    """

    if scene is None:
        scene = get_default_scene()

    if camera is None:
         camera = get_default_camera()

    if renderer is None:
        renderer = get_default_renderer(scene=scene, camera=camera)

    sources = parse_nodes(world)
    visualisers = parse_sources(sources)

    for _, visualiser in visualisers.items():

        visualiser.plot(renderer)

    if axes:
        renderer.scene.add(AxesHelper(0.1))

    return renderer, sources, visualisers