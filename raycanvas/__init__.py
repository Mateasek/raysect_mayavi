"""Ray tracing canvas for Cherab.

A visualization toolkit for Raysect scenegraphs using PyVista.
"""

__version__ = "0.1.0"
__author__ = "Cherab Development Team"
__license__ = "EUPL 1.1"

# Import main modules for easy access
from . import backend
from . import pyvista

__all__ = ["backend", "pyvista"]
