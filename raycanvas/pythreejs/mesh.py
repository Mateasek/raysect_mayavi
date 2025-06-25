from raycanvas.backend.builtin.primitives.mesh import TriangularMeshSource
from raycanvas.pythreejs.visualiser import PythreejsVisualiser

from pythreejs import (MeshStandardMaterial, MeshNormalMaterial, Group, BufferGeometry,
                       BufferAttribute, Mesh)


class TriangularMeshVisualiser(PythreejsVisualiser):
    """
    Class for visualising Raysect objetcs transformed to triangular meshses.
    """
    def __init__(self, source):

        
        self._init_plot_kwargs()

        
        super().__init__(source)

    def _init_object_kwargs(self):
        self.object_kwargs = {}
        self.object_kwargs["opacity"] = 0.5

    def set_source(self, source):

        if not isinstance(source, TriangularMeshSource):
            raise TypeError("source has to be instance of SourceBase.")
        
        self._source = source
    
    def get_data_object(self):

        # Define materials to be displayed
        material_wireframe = MeshStandardMaterial(color="#ffffff", wireframe=True)
        material_solid_back = MeshNormalMaterial(color="#888888",
                                            transparent=True,
                                            flatshading=True,
                                            opacity=0.5)
        material_solid_front = MeshNormalMaterial(color="#555555",
                                            transparent=True,
                                            flatshading=True,
                                            opacity=0.5)

        # Get vertices and triangles to build pythreejs mesh
        vertices = self._source.vertices

        # triangles have to be 1d vector of uint16
        triangles = self._source.triangles.astype("uint16").ravel()

        #group meshes
        group = Group()

        geometry = BufferGeometry(attributes=dict(
            position=BufferAttribute(vertices, normalized=False),
            index=BufferAttribute(triangles, normalized=False),
        #    color=BufferAttribute(vertexcolors),
        ))

        # create the visualisation of the front mesh faces
        front = Mesh(
            geometry=geometry,
            material=material_solid_front,
        )
        front.material.side = "FrontSide"
        front.renderorder = 0
        group.add(front)

        # create the visualisation of the back mesh faces
        back = Mesh(
            geometry=geometry,
            #material=MeshLambertMaterial(color='red', flatshading=False),
            material=material_solid_back,
        )
        back.material.side = "BackSide"
        back.renderorder = 1
        group.add(back)

        # create the visualisation of the mesh wires
        wires = Mesh(
            geometry=geometry,
            material=material_wireframe,
        )
        wires.renderorder = 2
        group.add(wires)

        return group

        
    def _add_object(self, plotter):
        
        mesh = self.get_data_object()
        return plotter.scene.add(mesh)