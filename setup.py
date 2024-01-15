from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import glob
import os

with open("VERSION", "r") as version_file:
        __version__ = version_file.read().strip()

build_flags = {}        
try:
    with open("src/Makevars", "r") as makevars:
        for line in makevars:
            parts = line.strip().split("=")
            if len(parts) == 2:
                build_flags[parts[0]] = parts[1]
except:
    build_flags = {"OPENMP_CXXFLAGS": "-fopenmp",
                   "OPENMP_LIBS": "-fopenmp",
                   "EIGEN_CXXFLAGS": "",
                   "EIGEN_LIBS": ""}
    
if os.getenv("CONDA_PREFIX") and build_flags["EIGEN_CXXFLAGS"] == "":
    ### conda
    build_flags["EIGEN_CXXFLAGS"] = "-I%s" % os.path.join(os.getenv("CONDA_PREFIX"), "include", "eigen3")
    
ext_modules = [
   Pybind11Extension(
        "scdemon_ext",
           sources=sorted(glob.glob("scdemon/py*.cpp")),
           include_dirs=["src"],
           cxx_std=17,
           define_macros=[("VERSION_INFO", __version__)],
           extra_compile_args=[build_flags["OPENMP_CXXFLAGS"], build_flags["EIGEN_CXXFLAGS"]],
           extra_link_args=[build_flags["OPENMP_LIBS"], build_flags["EIGEN_LIBS"]]
    ),
]

setup(
        name="scdemon",
        version=__version__,
        packages=find_packages(),
        ext_modules=ext_modules,
        cmdclass={"build_ext": build_ext}
)
