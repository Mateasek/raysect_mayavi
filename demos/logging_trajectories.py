
# External imports
# import matplotlib as mpl
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

# Raysect imports
from raysect.optical import World, translate, rotate, Point3D, d65_white, Ray, Vector3D
from raysect.optical.material.absorber import AbsorbingSurface
from raysect.optical.library import schott
from raysect.primitive import Sphere, Box
from raysect.optical.loggingray import LoggingRay
from raysect.primitive.lens.spherical import *


world = World()

# Create a glass BiConvex lens we want to study
lens_glass = schott("N-BK7")
lens_glass.transmission_only = True
lens = BiConvex(0.0254, 0.0052, 0.0506, 0.0506, parent=world, material=lens_glass)
lens.meta['viz-color'] = (66/255, 188/255, 244/255)
lens.meta['viz-opacity'] = 0.5

# lens.meta['viz']['color'] = (66/255, 188/255, 244/255)
# lens.meta['viz']['opacity'] = 0.5

# Create a target plane behind the lens.
target = Box(lower=Point3D(-0.05, -0.05, -0), upper=Point3D(0.05, 0.05, 0), material=AbsorbingSurface(),
             transform=translate(0, 0, 0.1), parent=world)
target.meta['viz-color'] = (224/255, 100/255, 17/255)


# for each sample direction trace a logging ray and plot the ray trajectory
plt.ion()
fig = plt.figure()
# ax = fig.gca(projection='3d')

# for u in np.linspace(-0.006, 0.006, 5):
for v in np.linspace(-0.012, 0.012, 11):

    start = Point3D(v, 0, -0.05)
    log_ray = LoggingRay(start, Vector3D(0, 0, 1))
    log_ray.trace(world)

    vertices = log_ray.path_vertices

    p = [(v.x, v.z) for v in vertices]
    p = np.array(p)

    plt.plot(p[:, 0], p[:, 1], 'k-')
    plt.plot(p[:, 0], p[:, 1], 'r.')


from raysect_mayavi import visualise_scenegraph
from mayavi import mlab


visualise_scenegraph(world)


for v in np.linspace(-0.012, 0.012, 11):

    start = Point3D(v, 0, -0.05)
    log_ray = LoggingRay(start, Vector3D(0, 0, 1))
    log_ray.trace(world)

    vertices = log_ray.path_vertices

    p = [(v.x, v.y, v.z) for v in vertices]
    p = np.array(p)

    mlab.plot3d(p[:, 0], p[:, 1], p[:, 2], tube_radius=0.0005)

plt.ioff()
plt.show()
