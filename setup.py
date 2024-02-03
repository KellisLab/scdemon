from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import glob
import os

with open("VERSION", "r") as version_file:
        __version__ = version_file.read().strip()

build_flags = {"OPENMP_CXXFLAGS": "-fopenmp",
               "OPENMP_LIBS": "-fopenmp",
               "EIGEN_CXXFLAGS": "",
               "EIGEN_LIBS": "",
               "GSL_CXXFLAGS": "",
               "GSL_LIBS": "-lgsl -lgslcblas -lcblas -lm".split()}
try:
    with open("src/Makevars", "r") as makevars:
        for line in makevars:
            parts = line.strip().split("=")
            if len(parts) == 2:
                build_flags[parts[0]] = parts[1]
except:
    pass

if os.getenv("CONDA_PREFIX"):
    ### conda
    build_flags["EIGEN_CXXFLAGS"] = "-I%s" % os.path.join(os.getenv("CONDA_PREFIX"), "include", "eigen3")

extra_compile_args = []
extra_link_args = []

def flag_compile(name):
    global extra_compile_args    
    x = build_flags.get(name, "")
    if isinstance(x, list):
        extra_compile_args += [y for y in x if y]
    elif x:
        extra_compile_args += [x]

def flag_link(name):
    global extra_link_args    
    x = build_flags.get(name, "")
    if isinstance(x, list):
        extra_link_args += [y for y in x if y]
    elif x:
        extra_link_args += [x]

flag_compile("OPENMP_CXXFLAGS")
flag_link("OPENMP_LIBS")
flag_compile("EIGEN_CXXFLAGS")
flag_link("EIGEN_LIBS")
flag_compile("GSL_CXXFLAGS")
flag_link("GSL_LIBS")
    
ext_modules = [
   Pybind11Extension(
        "scdemon_ext",
           sources=sorted(glob.glob("scdemon/py*.cpp")),
           include_dirs=["src"],
           cxx_std=17,
           define_macros=[("VERSION_INFO", __version__)],
           extra_compile_args=extra_compile_args,
           extra_link_args=extra_link_args
    ),
]

setup(
        name="scdemon",
        version=__version__,
        packages=find_packages(),
        ext_modules=ext_modules,
        cmdclass={"build_ext": build_ext}
)
