#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "se.hpp"

namespace py = pybind11;

PYBIND11_MODULE(scdemon_ext, m) {
    m.def("ols", &ols, "Ordinary least squares function");
}
