from collections import defaultdict
import sys
import os
import os.path as path
from pathlib import Path
import multiprocessing
import numpy
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

multiprocessing.set_start_method('fork')

force = False
profile = False
line_profile = False

if "--force" in sys.argv:
    force = True
    del sys.argv[sys.argv.index("--force")]

if "--profile" in sys.argv:
    profile = True
    del sys.argv[sys.argv.index("--profile")]

if "--line-profile" in sys.argv:
    line_profile = True
    del sys.argv[sys.argv.index("--line-profile")]

source_paths = ["raycanvas"]
compilation_includes = [".", numpy.get_include()]
compilation_args = ["-O3", "-Wno-unreachable-code-fallthrough"]
cython_directives = {"language_level": 3}
macros = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
setup_path = path.dirname(path.abspath(__file__))
num_processes = int(os.getenv("CHERAB_NCPU", "-1"))
if num_processes == -1:
    num_processes = multiprocessing.cpu_count()

if line_profile:
    compilation_args.append("-DCYTHON_TRACE=1")
    compilation_args.append("-DCYTHON_TRACE_NOGIL=1")
    cython_directives["linetrace"] = True
if profile:
    cython_directives["profile"] = True


extensions = []
for package in source_paths:
    for root, dirs, files in os.walk(path.join(setup_path, package)):
        for file in files:
            if path.splitext(file)[1] == ".pyx":
                pyx_file = path.relpath(path.join(root, file), setup_path)
                module = path.splitext(pyx_file)[0].replace("/", ".")
                extensions.append(
                    Extension(
                        module,
                        [pyx_file],
                        include_dirs=compilation_includes,
                        extra_compile_args=compilation_args,
                        define_macros=macros,
                    ),
                )


# generate .c files from .pyx
extensions = cythonize(
    extensions,
    nthreads=multiprocessing.cpu_count(),
    force=force,
    compiler_directives=cython_directives,
)

# Include demos in a separate directory in the distribution as data_files.
demo_parent_path = Path("share/cherab/demos/core")
data_files = defaultdict(list)
demos_source = Path("demos")
for item in demos_source.rglob("*"):
    if item.is_file():
        install_dir = demo_parent_path / item.parent.relative_to(demos_source)
        data_files[str(install_dir)].append(str(item))
data_files = list(data_files.items())


setup(
    name="raycanvas",
    version="0.1.0",
    license="EUPL 1.1",
    description="Ray tracing canvas for Cherab",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    install_requires=[
        "numpy>=1.21.0,<2.0",
        "scipy>=1.15.3",
        "pyvista[all,trame]",
        "jupyterlab",
        "raysect==0.8.1",
        "cherab>=1.5.0",
    ],
    extras_require={
        # Running ./dev/build_docs.sh runs setup.py, which requires cython.
        "docs": ["cython~=3.0", "sphinx", "sphinx-rtd-theme", "sphinx-tabs"],
        "test": ["pytest>=7.0.0", "pytest-cov>=4.0.0"],
    },
    packages=find_packages(include=["raycanvas*"]),
    package_data={"": [
        "**/*.pyx", "**/*.pxd",  # Needed to build Cython extensions.
        "**/*.json", "**/*.cl", "**/*.npy", "**/*.obj",  # Supplementary data
    ],
    },
    data_files=data_files,
    zip_safe=False,
    ext_modules=extensions,
    options=dict(
        build_ext={"parallel": num_processes},
    ),
)