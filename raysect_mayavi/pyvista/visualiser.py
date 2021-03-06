from raysect_mayavi.primitives.source import VisualiserBase

import pyvista as pv
pv.rcParams['use_ipyvtk'] = True


class PyvistaVisualiser(VisualiserBase):
    """
    This is the base class for the raysect_mayavi representation of raysect objects.
    """

    def __init__(self, source):

        self.set_source(source)
        self._init_figure_kwargs()
        self._init_plot_kwargs()

    def _init_figure_kwargs(self):
        self.figure_kwargs = {}
        self.figure_kwargs["window_size"] = (512, 256)

    def _init_plot_kwargs(self):
        self.plot_kwargs = {}

    @property
    def source(self):
        return self._source

    def set_source(self, source):
        raise NotImplementedError("Virtual method _set_source() has not been implemented.")

    def get_data_object(self):
        """
        Constructs the pyvista.DataObject representation of the Raysect primitive
        """
        raise NotImplementedError("Virtual method get_vista_object() has not been implemented.")

    def plot(self, plotter=None):
        """
        Plot the representation of the Raysect object into the pyvista.plotter.
        The representation is done always in the root node.

        :param plotter: Optional, specifies the pyvista.plotter to plot in.

        :return: pyvista.plotter
        """
        plotter = plotter or pv.Plotter(**self.figure_kwargs)
        if not isinstance(plotter, pv.Plotter):
            raise ValueError("plotter has to be instance of pyvista.Plotter.")

        self._add_object(plotter)

        return plotter

    def _add_object(self, plotter):
        raise NotImplementedError("Virtual method has not been implemented.")
