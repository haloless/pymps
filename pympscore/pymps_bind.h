#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using pyrvp = py::return_value_policy;

#define MPS_DEF_PYMOD(name, mod) void def_pymod_##name (py::module & mod)



