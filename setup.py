from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import glob

with open("VERSION", "r") as version_file:
        __version__ = version_file.read().strip()

ext_modules = [
   Pybind11Extension(
        "scdemon_ext",
           sources=sorted(glob.glob("scdemon/py*.cpp")),
           include_dirs=["src",
                         "/usr/include/eigen3",
                         "/usr/local/include/eigen3",
                         "../eigen",
                         ],
           cxx_std=17,
           define_macros=[("VERSION_INFO", __version__)],
           extra_compile_args=["-fopenmp"],
           extra_link_args=["-fopenmp"]
    ),
]

setup(
        name="scdemon",
        version=__version__,
        packages=find_packages(),
        ext_modules=ext_modules,
        cmdclass={"build_ext": build_ext}
)
