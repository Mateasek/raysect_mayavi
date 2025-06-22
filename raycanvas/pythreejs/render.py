from pythreejs import Scene, PerspectiveCamera, Renderer, OrbitControls, PointLight


def get_default_scene():
    return Scene()

def get_default_camera():
    return PerspectiveCamera(position=[1, 0, 0], fov=20)

def get_default_renderer(scene=None, camera=None):

    if scene is None:
        scene = get_default_scene()
    if camera is None:
        camera = get_default_camera()
        camera.children = [PointLight()]

    width = 800
    height = 600

    renderer = Renderer(camera=camera, scene=scene)
    renderer.controls = [OrbitControls(controlling=camera)]
    renderer.width = width
    renderer.height = height
    renderer.background = "black"
    renderer.background_opacity = 1

    return renderer





