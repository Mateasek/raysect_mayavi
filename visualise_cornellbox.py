{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "0faca6c8",
   "metadata": {},
   "outputs": [],
   "source": [
    "from raysect.optical import World\n",
    "from raysect.optical.loggingray import LoggingRay\n",
    "from raysect.core import Point3D, Vector3D, rotate_y\n",
    "\n",
    "from cornell_scene import add_cornell_box\n",
    "from raycanvas.pyvista import visualise"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "political-traffic",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "# prepare scene for visualisation\n",
    "world = World()\n",
    "camera, _ = add_cornell_box(world)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "functioning-olympus",
   "metadata": {},
   "outputs": [],
   "source": [
    "#visualise just the Cornell box\n",
    "plotter, visualisers = visualise(world)\n",
    "plotter.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7691ba76",
   "metadata": {},
   "outputs": [],
   "source": [
    "# add two logging rays to display ray paths\n",
    "origin = Point3D(0, 0, 0).transform(camera.to_root())\n",
    "direction = Vector3D(0, 0, 1).transform(camera.to_root())\n",
    "\n",
    "logray1 = LoggingRay(origin, direction)\n",
    "logray1.trace(world)\n",
    "\n",
    "logray2 = LoggingRay(origin, direction.transform(rotate_y(5)))\n",
    "logray2.trace(world)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "69180d98",
   "metadata": {},
   "outputs": [],
   "source": [
    "# add a logging ray to the plotter\n",
    "plotter, visualisers = visualise([logray1, logray2], plotter=plotter, visualisers=visualisers)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "52fc4b8e",
   "metadata": {},
   "outputs": [],
   "source": [
    "plotter.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "996ecb7d",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": ".venv",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.12.9"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
