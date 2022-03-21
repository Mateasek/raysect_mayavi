from raysect_mayavi.primitives.source import TriangularMeshSource
from raysect_mayavi.primitives.geometric.utility import triangulate_open_cylinder, triangulate_cylinder_lid, triangulate_spherical_cap
from raysect_mayavi.primitives.weld_vertices import weld_vertices

from raysect.primitive import Mesh
from raysect.primitive import Cylinder
from raysect.primitive.lens.spherical import BiConvex, BiConcave, PlanoConvex, PlanoConcave, Meniscus

import numpy as np
from math import sqrt


class CylindricalSource(TriangularMeshSource):
    """
    Class for creting mesh sources for objects with cylindrical symmetry e.g. cylinder and lenses. The mesh is composed
    of 3 surfaces: front surface, back surface and an open cylinder surface.
    """

    def __init__(self, raysect_object, vertical_divisions=10, cylindrical_divisions=36, radial_divisions=5):

        self._barrel_vertices = None
        self._barrel_triangles = None
        self._front_vertices = None
        self._front_triangles = None
        self._back_vertices = None
        self._back_triangles = None

        self.vertical_divisions = vertical_divisions
        self.cylindrical_divisions = cylindrical_divisions
        self.radial_divisions = radial_divisions

        super().__init__(raysect_object)

    def _graphic_source_from_raysect_object(self):

        # generate mesh parts if needed
        if self._back_triangles is None or self._back_vertices is None:
            self._generate_back_surface_mesh()
        if self._barrel_triangles is None or self._barrel_vertices is None:
            self._generate_barrel_surface_mesh()
        if self._front_triangles is None or self._front_vertices is None:
            self._generate_front_surface_mesh()

        #combine vertices of the 3 surfaces
        vertices = np.concatenate((self._back_vertices, self._barrel_vertices, self._front_vertices))

        # recalculate and combine triangle ids for merged vertices
        n_back = self._back_vertices.shape[0]
        n_barrel = self._barrel_vertices.shape[0]

        barrel_triangles = self._barrel_triangles.copy()
        barrel_triangles += n_back
        front_triangles = self._front_triangles.copy()
        front_triangles += n_back + n_barrel

        triangles = np.concatenate((self._back_triangles, barrel_triangles, front_triangles))

        # create raysect mesh and weld duplicit vertices
        mesh = Mesh(vertices=vertices, triangles=triangles)
        self._raysect_mesh = weld_vertices(mesh)

    def _generate_barrel_surface_mesh(self):
        """
        Generates triangular mesh for the barrel surface.
        :return:
        """
        raise NotImplementedError("Virtual method _generate_barrel_mesh() has not been implemented.")

    def _generate_front_surface_mesh(self):
        """
        Generates triangular mesh for the front surface.
        :return:
        """
        raise NotImplementedError("Virtual method _generate_front_surface_mesh() has not been implemented.")

    def _generate_back_surface_mesh(self):
        """
        Generates triangular mesh for the back surface.
        :return:
        """
        raise NotImplementedError("Virtual method _generate_back_surface_mesh() has not been implemented.")


class CylinderSource(CylindricalSource):
    """
    Class for graphical representation of the Raysect Cylinder primitive.
    :param raysect_object: Raysect Cone primitive instance
    :param vertical_divisions: Number of divisions of the barrel cylinder surface along the z axis (default is 10)
    :param cylindrical_divisions: Number of angular divisions of the cylinder surfaces (default is 36)
    :param radial_divisions: Nuber of radial divisions of the cylinder cap surfaces (default is 5)
    """

    def __init__(self, raysect_object, vertical_divisions=10, cylindrical_divisions=36, radial_divisions=5):

        if not isinstance(raysect_object, Cylinder):
            raise TypeError("The raysect_object has to be instance of Raysect Cylinder primitive, wrong type '{}' given.".format(type(raysect_object)))

        super().__init__(raysect_object, vertical_divisions, cylindrical_divisions, radial_divisions)

    def _generate_barrel_surface_mesh(self):

        radius = self.raysect_object.radius
        height = self.raysect_object.height

        vertices, triangles = triangulate_open_cylinder(radius, height, self.cylindrical_divisions,
                                                        self.vertical_divisions)

        self._barrel_vertices = vertices
        self._barrel_triangles = triangles

    def _generate_front_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = self.raysect_object.radius
        edge = self._barrel_vertices[0:self.cylindrical_divisions, 0:2]

        vertices, triangles = triangulate_cylinder_lid(radius, self.radial_divisions, edge_vertices=edge)

        vertices[:, 2] = self.raysect_object.height

        self._front_vertices = vertices
        self._front_triangles = triangles

    def _generate_back_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = self.raysect_object.radius
        edge = self._barrel_vertices[0:self.cylindrical_divisions, 0:2]

        vertices, triangles = triangulate_cylinder_lid(radius, self.radial_divisions, edge_vertices=edge)

        self._back_vertices = vertices
        self._back_triangles = triangles


class BiConvexLensSource(CylindricalSource):
    """
    Class for graphical representation of the Raysect BiConvex lens primitive.
    :param raysect_object: Raysect Biconvex primitive instance
    :param vertical_divisions: Number of divisions of the barrel cylinder surface along the z axis (default is 10)
    :param cylindrical_divisions: Number of angular divisions of the lens surfaces (default is 36)
    :param radial_divisions: Nuber of radial divisions of the front and back lens surfaces (default is 5)
    """

    def __init__(self, raysect_object, vertical_divisions=10, cylindrical_divisions=36, radial_divisions=5):

        if not isinstance(raysect_object, BiConvex):
            raise TypeError("The raysect_object has to be instance of Raysect BiConvex primitive, wrong type '{}' given.".format(type(raysect_object)))

        super().__init__(raysect_object, vertical_divisions, cylindrical_divisions, radial_divisions)

    def _generate_barrel_surface_mesh(self):

        radius = 0.5 * self.raysect_object.diameter
        height = self.raysect_object.edge_thickness

        vertices, triangles = triangulate_open_cylinder(radius, height, self.cylindrical_divisions,
                                                        self.vertical_divisions)

        # shift cylinder vertices z coordinate
        vertices[:, 2] += self.raysect_object.back_thickness

        self._barrel_vertices = vertices
        self._barrel_triangles = triangles

    def _generate_front_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[-self.cylindrical_divisions::, 0:2]
        z_center = self.raysect_object.center_thickness - self.raysect_object.front_curvature
        curvature = self.raysect_object.front_curvature

        vertices, triangles = triangulate_spherical_cap(curvature, radius,
                                                        self.radial_divisions, edge_vertices=edge)

        for i, v in enumerate(vertices):
            r2 = v[0] ** 2 + v[1] ** 2
            z = z_center + sqrt(curvature ** 2 - r2)
            vertices[i, 2] = z

        self._front_vertices = vertices
        self._front_triangles = triangles

    def _generate_back_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[0:self.cylindrical_divisions, 0:2]
        z_center = self.raysect_object.back_curvature
        curvature = self.raysect_object.back_curvature

        vertices, triangles = triangulate_spherical_cap(curvature, radius, self.radial_divisions, edge_vertices=edge)

        for i, v in enumerate(vertices):
            r2 = v[0] ** 2 + v[1] ** 2
            z = z_center - sqrt(curvature ** 2 - r2)
            vertices[i, 2] = z

        self._back_vertices = vertices
        self._back_triangles = triangles


class BiConcaveLensSource(CylindricalSource):
    """
    Class for graphical representation of the Raysect BiConcave lens primitive.
    :param raysect_object: Raysect BiConcave primitive instance
    :param vertical_divisions: Number of divisions of the barrel cylinder surface along the z axis (default is 10)
    :param cylindrical_divisions: Number of angular divisions of the lens surfaces (default is 36)
    :param radial_divisions: Nuber of radial divisions of the front and back lens surfaces (default is 5)
    """

    def __init__(self, raysect_object, vertical_divisions=10, cylindrical_divisions=36, radial_divisions=5):

        if not isinstance(raysect_object, BiConcave):
            raise TypeError("The raysect_object has to be instance of Raysect BiConcave primitive, wrong type '{}' given.".format(type(raysect_object)))

        super().__init__(raysect_object, vertical_divisions, cylindrical_divisions, radial_divisions)

    def _generate_barrel_surface_mesh(self):

        radius = 0.5 * self.raysect_object.diameter
        height = self.raysect_object.edge_thickness

        vertices, triangles = triangulate_open_cylinder(radius, height, self.cylindrical_divisions,
                                                        self.vertical_divisions)

        # shift cylinder vertices z coordinate
        vertices[:, 2] -= self.raysect_object.back_thickness

        self._barrel_vertices = vertices
        self._barrel_triangles = triangles

    def _generate_front_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[-self.cylindrical_divisions::, 0:2]
        z_center = self.raysect_object.center_thickness + self.raysect_object.front_curvature
        curvature = self.raysect_object.front_curvature

        vertices, triangles = triangulate_spherical_cap(curvature, radius,
                                                        self.radial_divisions, edge_vertices=edge)

        for i, v in enumerate(vertices):
            r2 = v[0] ** 2 + v[1] ** 2
            z = z_center - sqrt(curvature ** 2 - r2)
            vertices[i, 2] = z

        self._front_vertices = vertices
        self._front_triangles = triangles

    def _generate_back_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[0:self.cylindrical_divisions, 0:2]
        z_center = -self.raysect_object.back_curvature
        curvature = self.raysect_object.back_curvature

        vertices, triangles = triangulate_spherical_cap(curvature, radius, self.radial_divisions, edge_vertices=edge)

        for i, v in enumerate(vertices):
            r2 = v[0] ** 2 + v[1] ** 2
            z = z_center + sqrt(curvature ** 2 - r2)
            vertices[i, 2] = z

        self._back_vertices = vertices
        self._back_triangles = triangles


class PlanoConvexLensSource(CylindricalSource):
    """
    Class for graphical representation of the Raysect PlanoConvex lens primitive.
    :param raysect_object: Raysect PlanoConvex primitive instance
    :param vertical_divisions: Number of divisions of the barrel cylinder surface along the z axis (default is 10)
    :param cylindrical_divisions: Number of angular divisions of the lens surfaces (default is 36)
    :param radial_divisions: Nuber of radial divisions of the front and back lens surfaces (default is 5)
    """

    def __init__(self, raysect_object, vertical_divisions=10, cylindrical_divisions=36, radial_divisions=5):

        if not isinstance(raysect_object, PlanoConvex):
            raise TypeError("The raysect_object has to be instance of Raysect PlanoConvex primitive, wrong type '{}' given.".format(type(raysect_object)))

        super().__init__(raysect_object, vertical_divisions, cylindrical_divisions, radial_divisions)

    def _generate_barrel_surface_mesh(self):

        radius = 0.5 * self.raysect_object.diameter
        height = self.raysect_object.edge_thickness

        vertices, triangles = triangulate_open_cylinder(radius, height, self.cylindrical_divisions,
                                                        self.vertical_divisions)

        self._barrel_vertices = vertices
        self._barrel_triangles = triangles

    def _generate_front_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[-self.cylindrical_divisions::, 0:2]
        z_center = self.raysect_object.center_thickness - self.raysect_object.curvature
        curvature = self.raysect_object.curvature

        vertices, triangles = triangulate_spherical_cap(curvature, radius,
                                                        self.radial_divisions, edge_vertices=edge)

        for i, v in enumerate(vertices):
            r2 = v[0] ** 2 + v[1] ** 2
            z = z_center + sqrt(curvature ** 2 - r2)
            vertices[i, 2] = z

        self._front_vertices = vertices
        self._front_triangles = triangles

    def _generate_back_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[0:self.cylindrical_divisions, 0:2]

        vertices, triangles = triangulate_cylinder_lid(radius, self.radial_divisions, edge_vertices=edge)

        self._back_vertices = vertices
        self._back_triangles = triangles


class PlanoConcaveLensSource(CylindricalSource):
    """
    Class for graphical representation of the Raysect PlanoConcave lens primitive.
    :param raysect_object: Raysect PlanoConcave primitive instance
    :param vertical_divisions: Number of divisions of the barrel cylinder surface along the z axis (default is 10)
    :param cylindrical_divisions: Number of angular divisions of the lens surfaces (default is 36)
    :param radial_divisions: Nuber of radial divisions of the front and back lens surfaces (default is 5)
    """

    def __init__(self, raysect_object, vertical_divisions=10, cylindrical_divisions=36, radial_divisions=5):

        if not isinstance(raysect_object, PlanoConcave):
            raise TypeError("The raysect_object has to be instance of Raysect PlanoConcave primitive, wrong type '{}' given.".format(type(raysect_object)))

        super().__init__(raysect_object, vertical_divisions, cylindrical_divisions, radial_divisions)

    def _generate_barrel_surface_mesh(self):

        radius = 0.5 * self.raysect_object.diameter
        height = self.raysect_object.edge_thickness

        vertices, triangles = triangulate_open_cylinder(radius, height, self.cylindrical_divisions,
                                                        self.vertical_divisions)

        self._barrel_vertices = vertices
        self._barrel_triangles = triangles

    def _generate_front_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[-self.cylindrical_divisions::, 0:2]
        z_center = self.raysect_object.center_thickness + self.raysect_object.curvature
        curvature = self.raysect_object.curvature

        vertices, triangles = triangulate_spherical_cap(curvature, radius,
                                                        self.radial_divisions, edge_vertices=edge)

        for i, v in enumerate(vertices):
            r2 = v[0] ** 2 + v[1] ** 2
            z = z_center - sqrt(curvature ** 2 - r2)
            vertices[i, 2] = z

        self._front_vertices = vertices
        self._front_triangles = triangles

    def _generate_back_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[0:self.cylindrical_divisions, 0:2]

        vertices, triangles = triangulate_cylinder_lid(radius, self.radial_divisions, edge_vertices=edge)

        self._back_vertices = vertices
        self._back_triangles = triangles


class MeniscusLensSource(CylindricalSource):
    """
    Class for graphical representation of the Raysect Meniscus lens primitive.
    :param raysect_object: Raysect Meniscus primitive instance
    :param vertical_divisions: Number of divisions of the barrel cylinder surface along the z axis (default is 10)
    :param cylindrical_divisions: Number of angular divisions of the lens surfaces (default is 36)
    :param radial_divisions: Nuber of radial divisions of the front and back lens surfaces (default is 5)
    """

    def __init__(self, raysect_object, vertical_divisions=10, cylindrical_divisions=36, radial_divisions=5):

        if not isinstance(raysect_object, Meniscus):
            raise TypeError("The raysect_object has to be instance of Raysect Meniscus primitive, wrong type '{}' given.".format(type(raysect_object)))

        super().__init__(raysect_object, vertical_divisions, cylindrical_divisions, radial_divisions)

    def _generate_barrel_mesh(self):

        radius = 0.5 * self.raysect_object.diameter
        height = self.raysect_object.edge_thickness

        vertices, triangles = triangulate_open_cylinder(radius, height, self.cylindrical_divisions,
                                                        self.vertical_divisions)

        # shift cylinder vertices z coordinate
        vertices[:, 2] -= self.raysect_object.back_thickness

        self._barrel_vertices = vertices
        self._barrel_triangles = triangles

    def _generate_front_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[-self.cylindrical_divisions::, 0:2]
        z_center = self.raysect_object.center_thickness - self.raysect_object.front_curvature
        curvature = self.raysect_object.front_curvature

        vertices, triangles = triangulate_spherical_cap(curvature, radius,
                                                        self.radial_divisions, edge_vertices=edge)

        for i, v in enumerate(vertices):
            r2 = v[0] ** 2 + v[1] ** 2
            z = z_center + sqrt(curvature ** 2 - r2)
            vertices[i, 2] = z

        self._front_vertices = vertices
        self._front_triangles = triangles

    def _generate_back_surface_mesh(self):

        if self._barrel_vertices is None or self._barrel_triangles is None:
            self._generate_barrel_surface_mesh()

        radius = 0.5 * self.raysect_object.diameter
        edge = self._barrel_vertices[0:self.cylindrical_divisions, 0:2]
        z_center = -self.raysect_object.back_curvature
        curvature = self.raysect_object.back_curvature

        vertices, triangles = triangulate_spherical_cap(curvature, radius, self.radial_divisions, edge_vertices=edge)

        for i, v in enumerate(vertices):
            r2 = v[0] ** 2 + v[1] ** 2
            z = z_center + sqrt(curvature ** 2 - r2)
            vertices[i, 2] = z

        self._back_vertices = vertices
        self._back_triangles = triangles

    def _generate_barrel_surface_mesh(self):

        radius = 0.5 * self.raysect_object.diameter
        height = self.raysect_object.edge_thickness

        vertices, triangles = triangulate_open_cylinder(radius, height, self.cylindrical_divisions,
                                                        self.vertical_divisions)

        # shift cylinder vertices z coordinate
        vertices[:, 2] -= self.raysect_object.back_thickness

        self._barrel_vertices = vertices
        self._barrel_triangles = triangles