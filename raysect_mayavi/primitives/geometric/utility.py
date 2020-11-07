import numpy as np

from math import sqrt

from scipy.spatial import Delaunay


def triangulate_circle(radius, radial_divisions=10, edge_vertices=None, cylindrical_divisions=90):
    """
    Meshes a cirlcle using scipy.spatial.Delaunay. Mesh vertices are generated on co-centric circles.
    :param radius: Radius of the circle to mesh
    :param radial_divisions: Sets the number of co-centric sub-circles the circle is separated into and
                              vertices are placed at.
    :param edge_vertices: 2D coordinates (x, y) of vertices on the outer-most circle to be used. This is useful if the
                           generated mesh is supposed to be combined with some other mesh parts to generate a more
                           complicated mesh object.
    :param cylindrical_divisions: Number of vertices the outer-most circle will be divided in. Is ignored if
                                   outer_vertices is passed.
    :return: tupple (vertices, triangles) where vertices and triangles are 2D numpy arrays with vertex coordinates and
             triangle vertex ids, respectively.
    """
    if edge_vertices is not None:
        vertices = edge_vertices.tolist()
        cylindrical_divisions = len(edge_vertices)
        working_cylindrical_segments = cylindrical_divisions - int(cylindrical_divisions / (radial_divisions + 1))
        segment_start = 1
    else:
        working_cylindrical_segments = cylindrical_divisions
        vertices = []
        segment_start = 0

    for i in range(segment_start, radial_divisions):

        working_radius = radius * (1 - (i / radial_divisions))
        theta_step = 360 / working_cylindrical_segments

        for j in range(working_cylindrical_segments):
            theta_rad = np.deg2rad(j * theta_step)
            vertices.append([working_radius * np.cos(theta_rad), working_radius * np.sin(theta_rad)])

        working_cylindrical_segments -= int(cylindrical_divisions / (radial_divisions + 1))
        if working_cylindrical_segments < 5:
            working_cylindrical_segments = 5
    # Finally, add centre point and convert to numpy array
    vertices.append([0, 0])
    vertices = np.array(vertices)

    triangles = Delaunay(vertices).simplices

    return vertices, triangles


def triangulate_open_cylinder(radius, height, cylindrical_divisions=90, vertical_divisions=10):
    """
    Generates a 3D mesh of an opene cylinder. The cylindrer axis as aligned with the z axis and extends from z=0 to
    z=height.
    :param cylindrical_divisions: Number of cylindrical division of the barrel the vertices will be generated on.
    :param vertical_divisions: Number of vertical levels vertices will be generated on.
    :param radius: Radius of the barrel.
    :param height: Height of the barrel.
    :return: tupple (vertices, triangles) where vertices is a 2D numpy array with 3D (Mx3) vertex coorditanes. The
             triangles is a 2D numpy array with triangle vertex ids.
    """
    theta_step = 360 / cylindrical_divisions

    vertices = []
    for i in range(vertical_divisions):
        z = (i / (vertical_divisions - 1)) * height
        for j in range(cylindrical_divisions):
            theta_rad = np.deg2rad(j * theta_step)
            vertices.append([radius * np.cos(theta_rad), radius * np.sin(theta_rad), z])

    triangles = []
    for i in range(vertical_divisions - 1):

        row_start = cylindrical_divisions * i
        next_row_start = cylindrical_divisions * (i + 1)

        for j in range(cylindrical_divisions):

            v1 = row_start + j
            if j != cylindrical_divisions - 1:
                v2 = row_start + j + 1
                v3 = next_row_start + j + 1
            else:
                v2 = row_start
                v3 = next_row_start
            v4 = next_row_start + j

            triangles.append([v1, v2, v3])
            triangles.append([v3, v4, v1])

    # convert to numpy arrays
    vertices = np.array(vertices)
    triangles = np.array(triangles)

    return vertices, triangles


def triangulate_cylinder_lid(radius, radial_divisions, edge_vertices=None, cylindrical_divisions=90):
    """

    :param radial_divisions:
    :param radius:
    :param z_height:
    :param edge_vertices:
    :return:
    """

    circle_vertices, triangles =triangulate_circle(radius, radial_divisions, edge_vertices=edge_vertices,
                                                   cylindrical_divisions=cylindrical_divisions)

    vertices = np.zeros((circle_vertices.shape[0], 3))
    vertices[:, 0:-1] = circle_vertices

    return vertices, triangles


def triangulate_spherical_cap(sphere_radius, base_radius, radial_divisions, edge_vertices=None, cylindrical_divisions=90):
    """
    Triangulate a spherical cap without the base with axis of symmetry aligned with the z axis and the sphere center at
    z=0. The cap is located in the positive z half-plane.
    :param sphere_radius: Radius of the sphere
    :param base_radius: The radius of the cap base
    :param radial_divisions: Number of radial segments to place vertices at
    :param edge_vertices: Gives the edge vertex positions (x, y). This is useful if the cap edge points should be identical to
                          other meshes which simplifies welding of meshes.
    :return: tupple (vertices, triangles) where vertices is a 2D numpy array with 3D (Mx3) vertex coorditanes. The
             triangles is a 2D numpy array with triangle vertex ids.
    """""

    circle_vertices, triangles = triangulate_circle(base_radius, radial_divisions, edge_vertices=edge_vertices,
                                                    cylindrical_divisions=cylindrical_divisions)

    # lift circle vertices to form the 3D surface
    vertices = np.zeros((circle_vertices.shape[0], 3))
    vertices[:, 0:-1] = circle_vertices
    sr2 = sphere_radius ** 2
    for i, v in enumerate(vertices):
        r2 = v[0] ** 2 + v[1] ** 2
        z = sqrt(sr2 - r2)
        vertices[i, -1] = z


    return vertices, triangles

def triangulate_parabolic_cap(parabola_height, base_radius, radial_divisions, edge_vertices=None, cylindrical_divisions=90):
    """
    The parabolic cap is defined by a radius and height. It lies along the z-axis and extends over the z range [0, height].
     The base of the parabola is opened. The base of the parabola lies on the x-y plane, the parabola vertex (tip) lies 
     at z=height.    
    :param parabola_height: Intersection of the parabola with the z axis.
    :param base_radius: The radius of the cap base
    :param radial_divisions: Number of radial segments to place vertices at
    :param edge_vertices: Gives the edge vertex positions (x, y). This is useful if the cap edge points should be identical to
                          other meshes which simplifies welding of meshes.
    :return: tupple (vertices, triangles) where vertices is a 2D numpy array with 3D (Mx3) vertex coorditanes. The
             triangles is a 2D numpy array with triangle vertex ids.
    """""

    circle_vertices, triangles = triangulate_circle(base_radius, radial_divisions, edge_vertices=edge_vertices,
                                                    cylindrical_divisions=cylindrical_divisions)

    br2 = base_radius ** 2

    # lift circle vertices to form the 3D surface
    vertices = np.zeros((circle_vertices.shape[0], 3))
    vertices[:, 0:-1] = circle_vertices
    for i, v in enumerate(vertices):
        r2 = v[0] ** 2 + v[1] ** 2
        z = parabola_height * (1 - r2 / br2)
        vertices[i, -1] = z

    return vertices, triangles

def triangulate_conical_cap(cone_height, base_radius, radial_divisions, edge_vertices=None, cylindrical_divisions=90):
    """
    The conical cap is defined by a radius and height. It lies along the z-axis and extends over the z range [0, height].
     The base of the cylinder is opened. The base of the cylinder lies on the x-y plane, the cone vertex (tip) lies 
     at z=height.    
    :param cone_height: Height of the cone, intersection of the conical cap with the z axis.
    :param base_radius: The radius of the cap base
    :param radial_divisions: Number of radial segments to place vertices at
    :param edge_vertices: Gives the edge vertex positions (x, y). This is useful if the cap edge points should be identical to
                          other meshes which simplifies welding of meshes.
    :return: tupple (vertices, triangles) where vertices is a 2D numpy array with 3D (Mx3) vertex coorditanes. The
             triangles is a 2D numpy array with triangle vertex ids.
    """""

    circle_vertices, triangles = triangulate_circle(base_radius, radial_divisions, edge_vertices=edge_vertices,
                                                    cylindrical_divisions=cylindrical_divisions)
    
    z_height = lambda x: cone_height * (1 - x / base_radius)
    
    # lift circle vertices to form the 3D surface
    vertices = np.zeros((circle_vertices.shape[0], 3))
    vertices[:, 0:-1] = circle_vertices
    
    for i, v in enumerate(vertices):
        r = sqrt(v[0] ** 2 + v[1] ** 2)
        z = z_height(r)
        vertices[i, -1] = z
    
    return vertices, triangles
        
        