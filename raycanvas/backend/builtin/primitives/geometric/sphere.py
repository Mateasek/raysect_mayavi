import numpy as np

from raycanvas.backend.base.source import TriangularMeshSource

from raysect.core import Vector3D
from raysect.primitive import Sphere, Mesh

class SphereSource(TriangularMeshSource):
    """
    Class for graphical representation of the Raysect Sphere primitive.
    :param raysect_object: Raysect Sphere primitive instance
    :param subdivision_count: Number of subdivisions to perform on the basic icosohedren with 12 vertices and 20 faces
    """

    def __init__(self, raysect_object, subdivision_count=2):

        if not isinstance(raysect_object, Sphere):
            raise TypeError("The raysect_object has to be instance of Raysect Sphere primitive, wrong type '{}' given.".format(type(raysect_object)))

        self.subdivision_count = subdivision_count

        super().__init__(raysect_object)

    def _graphic_source_from_raysect_object(self):

        # Calculate vertices and faces using the icosohedren method
        # We compute a regular icosohedren with 12 vertices and 20 faces.
        # Vertices given by all perturbations of:
        # (0, ±1, ±ϕ), (±1, ±ϕ, 0), (±ϕ, 0, ±1), where ϕ = golden ratio

        golden_ratio = 1.61803398875

        radius = self._raysect_object.radius

        v1 = Vector3D(-1.0, golden_ratio, 0.0).normalise() * radius
        v2 = Vector3D(1.0, golden_ratio, 0.0).normalise() * radius
        v3 = Vector3D(-1.0, -golden_ratio, 0.0).normalise() * radius
        v4 = Vector3D(1.0, -golden_ratio, 0.0).normalise() * radius
        v5 = Vector3D(0.0, -1.0, golden_ratio).normalise() * radius
        v6 = Vector3D(0.0, 1.0, golden_ratio).normalise() * radius
        v7 = Vector3D(0.0, -1.0, -golden_ratio).normalise() * radius
        v8 = Vector3D(0.0, 1.0, -golden_ratio).normalise() * radius
        v9 = Vector3D(golden_ratio, 0.0, -1.0).normalise() * radius
        v10 = Vector3D(golden_ratio, 0.0, 1.0).normalise() * radius
        v11 = Vector3D(-golden_ratio, 0.0, -1.0).normalise() * radius
        v12 = Vector3D(-golden_ratio, 0.0, 1.0).normalise() * radius

        vertices = [
            [v1.x, v1.y, v1.z],
            [v2.x, v2.y, v2.z],
            [v3.x, v3.y, v3.z],
            [v4.x, v4.y, v4.z],
            [v5.x, v5.y, v5.z],
            [v6.x, v6.y, v6.z],
            [v7.x, v7.y, v7.z],
            [v8.x, v8.y, v8.z],
            [v9.x, v9.y, v9.z],
            [v10.x, v10.y, v10.z],
            [v11.x, v11.y, v11.z],
            [v12.x, v12.y, v12.z],
        ]

        triangles = [
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1]
        ]

        # Optional - subdivision of icosohedren to increase resolution
        num_vertices = 12
        num_triangles = 20
        for i in range(self._subdivision_count):
            for j in range(num_triangles):
                triangle = triangles[j]
                # extract current triangle vertices
                v0_id = triangle[0]
                v1_id = triangle[1]
                v2_id = triangle[2]
                v0 = Vector3D(vertices[v0_id][0], vertices[v0_id][1], vertices[v0_id][2])
                v1 = Vector3D(vertices[v1_id][0], vertices[v1_id][1], vertices[v1_id][2])
                v2 = Vector3D(vertices[v2_id][0], vertices[v2_id][1], vertices[v2_id][2])

                # subdivide with three new vertices
                v3 = (v0 + v1).normalise() * radius
                v3_id = num_vertices
                v4 = (v1 + v2).normalise() * radius
                v4_id = num_vertices + 1
                v5 = (v2 + v0).normalise() * radius
                v5_id = num_vertices + 2
                vertices.append([v3.x, v3.y, v3.z])
                vertices.append([v4.x, v4.y, v4.z])
                vertices.append([v5.x, v5.y, v5.z])

                # ... and three new faces
                triangles[j] = [v0_id, v3_id, v5_id]  # replace the first face
                triangles.append([v3_id, v1_id, v4_id])
                triangles.append([v4_id, v2_id, v5_id])
                triangles.append([v3_id, v4_id, v5_id])

                num_vertices += 3
                num_triangles += 3

        vertices = np.array(vertices)
        triangles = np.array(triangles)

        self._raysect_mesh = Mesh(vertices, triangles)

    @property
    def subdivision_count(self):
        return self._subdivision_count

    @subdivision_count.setter
    def subdivision_count(self, value):

        value = int(value)

        if not value > 0:
            raise ValueError("value has to be larger than 0")

        self._subdivision_count = value
