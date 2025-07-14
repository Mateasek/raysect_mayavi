from raysect.primitive import Mesh, Box, Sphere, Cylinder, Cone, Parabola, Intersect, Union, Subtract
from raysect.primitive.csg import CSGPrimitive
from raysect.primitive.lens.spherical import BiConvex, BiConcave, PlanoConvex, PlanoConcave, Meniscus
from raysect.optical.loggingray import LoggingRay
from raysect.optical.observer import SightLine, FibreOptic, Pixel, TargettedPixel, MeshPixel, MeshCamera
from raysect.optical.observer import CCDArray, TargettedCCDArray, OrthographicCamera, PinholeCamera, VectorCamera

from cherab.tools.observers.spectroscopy import SpectroscopicFibreOptic

from raycanvas.backend.builtin.primitives.mesh_tools import subdivide
from raycanvas.backend.builtin.primitives.mesh_csg import perform_mesh_csg
from raycanvas.backend.builtin.primitives.mesh_csg import (
    Intersect as IntersectOperator,
    Union as UnionOperator,
    Subtract as SubtractOperator,
)
from raycanvas.backend.builtin.primitives.geometric import BoxSource, SphereSource, CylinderSource, ConeSource, ParabolaSource
from raycanvas.backend.builtin.primitives.geometric import (
    BiConvexLensSource,
    BiConcaveLensSource,
    PlanoConvexLensSource,
    PlanoConcaveLensSource,
    MeniscusLensSource,
)
from raycanvas.backend.base.ray import LoggingRaySource
from raycanvas.backend.base.source import TriangularMeshSource
from raycanvas.backend.base.observer import ObserverSource
from raycanvas.backend.builtin.primitives.mesh_tools import subdivide


class MeshSource(TriangularMeshSource):
    """
    Class for graphical representation of the Raysect Mesh primitive.
    :param raysect_object: Raysect Mesh primitive instance
    """

    def __init__(self, raysect_object):
        if not isinstance(raysect_object, Mesh):
            raise TypeError(
                "The raysect_object has to be instance of Raysect " "Box primitive, wrong type '{}' given.".format(
                    type(raysect_object)
                )
            )

        super().__init__(raysect_object)

    def _graphic_source_from_raysect_object(self):
        self._raysect_mesh = self._raysect_object


class CSGMeshSource(TriangularMeshSource):
    """
    Class for graphical representation of the Raysect CSGPrimitive primitive.
    :param raysect_object: Raysect CSGPrimitive primitive instance
    """

    def __init__(self, raysect_object):
        if not isinstance(raysect_object, CSGPrimitive):
            raise TypeError(
                "The raysect_object has to be instance of Raysect Box primitive, wrong type '{}' given.".format(
                    type(raysect_object)
                )
            )

        super().__init__(raysect_object)

    def _graphic_source_from_raysect_object(self):
        primitive_a = self._raysect_object.primitive_a
        primitive_b = self._raysect_object.primitive_b

        mesh_a = to_mesh(primitive_a)
        mesh_a_transformed = Mesh(mesh_a._vertices_transform(primitive_a.transform), mesh_a.triangles)
        mesh_b = to_mesh(primitive_b)
        mesh_b_transformed = Mesh(mesh_b._vertices_transform(primitive_b.transform), mesh_b.triangles)

        if self._raysect_object.__class__ == Intersect:
            self._csg_operator = IntersectOperator()
        elif self._raysect_object.__class__ == Union:
            self._csg_operator = UnionOperator()
        elif self._raysect_object.__class__ == Subtract:
            self._csg_operator = SubtractOperator()
        else:
            raise ValueError("Unidentified CSG primitive '{}'.".format(self._raysect_object.__class__))

        self._raysect_mesh = perform_mesh_csg(mesh_a_transformed, mesh_b_transformed, operator=self._csg_operator)


_object_handlers = {
    Box: BoxSource,
    Sphere: SphereSource,
    Cylinder: CylinderSource,
    Cone: ConeSource,
    Mesh: MeshSource,
    Intersect: CSGMeshSource,
    Union: CSGMeshSource,
    Subtract: CSGMeshSource,
    BiConvex: BiConvexLensSource,
    BiConcave: BiConcaveLensSource,
    PlanoConvex: PlanoConvexLensSource,
    PlanoConcave: PlanoConcaveLensSource,
    Meniscus: MeniscusLensSource,
    Parabola: ParabolaSource,
    LoggingRay: LoggingRaySource,
    SightLine: ObserverSource,
    FibreOptic: ObserverSource,
    Pixel: ObserverSource,
    TargettedPixel: ObserverSource,
    MeshPixel: ObserverSource,
    MeshCamera: ObserverSource,
    CCDArray: ObserverSource,
    TargettedCCDArray: ObserverSource,
    OrthographicCamera: ObserverSource,
    PinholeCamera: ObserverSource,
    VectorCamera: ObserverSource,
    SpectroscopicFibreOptic: ObserverSource,
}


def to_mesh(primitive):
    """
    Automatically assings graphical representation to the primitive. Primitive has to be any of Raysect
    primitives or observers.

    :param primitive: Raysect primitive or observer object.
    :return MayaviSource
    """
    for raysect_class in _object_handlers.keys(): 
        if isinstance(primitive, raysect_class):
            handler = _object_handlers[raysect_class]
            return handler(primitive)

    raise ValueError("Unrecognised Raysect primitive, '{}'.".format(type(primitive)))
