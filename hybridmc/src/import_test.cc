#include <pybind11/embed.h>
#include <iostream>

namespace py = pybind11;

int main() {
    py::scoped_interpreter python;

    py::module_ sklearn = py::module_::import("sklearn");
    py::print(sklearn.attr("linear_model"));

}

