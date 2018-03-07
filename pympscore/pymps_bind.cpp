#if 0
#include <Python.h>

/*
 * Implements an example function.
 */
PyDoc_STRVAR(pympscore_example_doc, "example(obj, number)\
\
Example function");

PyObject *pympscore_example(PyObject *self, PyObject *args, PyObject *kwargs) {
    /* Shared references that do not need Py_DECREF before returning. */
    PyObject *obj = NULL;
    int number = 0;

    /* Parse positional and keyword arguments */
    static char* keywords[] = { "obj", "number", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", keywords, &obj, &number)) {
        return NULL;
    }

    /* Function implementation starts here */

    if (number < 0) {
        PyErr_SetObject(PyExc_ValueError, obj);
        return NULL;    /* return NULL indicates error */
    }

    Py_RETURN_NONE;
}

/*
 * List of functions to add to pympscore in exec_pympscore().
 */
static PyMethodDef pympscore_functions[] = {
    { "example", (PyCFunction)pympscore_example, METH_VARARGS | METH_KEYWORDS, pympscore_example_doc },
    { NULL, NULL, 0, NULL } /* marks end of array */
};

/*
 * Initialize pympscore. May be called multiple times, so avoid
 * using static state.
 */
int exec_pympscore(PyObject *module) {
    PyModule_AddFunctions(module, pympscore_functions);

    PyModule_AddStringConstant(module, "__author__", "sun");
    PyModule_AddStringConstant(module, "__version__", "1.0.0");
    PyModule_AddIntConstant(module, "year", 2018);

    return 0; /* success */
}

/*
 * Documentation for pympscore.
 */
PyDoc_STRVAR(pympscore_doc, "The pympscore module");


static PyModuleDef_Slot pympscore_slots[] = {
    { Py_mod_exec, exec_pympscore },
    { 0, NULL }
};

static PyModuleDef pympscore_def = {
    PyModuleDef_HEAD_INIT,
    "pympscore",
    pympscore_doc,
    0,              /* m_size */
    NULL,           /* m_methods */
    pympscore_slots,
    NULL,           /* m_traverse */
    NULL,           /* m_clear */
    NULL,           /* m_free */
};

PyMODINIT_FUNC PyInit_pympscore() {
    return PyModuleDef_Init(&pympscore_def);
}
#endif

////////////////////////////////////////////////////////////////////////////////

#include "pymps_bind.h"

#pragma warning(push)
#pragma warning(disable: 4996)

////////////////////////////////////////////////////////////////////////////////


//--------------------------------------------------------------------------------

//
// separate binding codes to reduce compile time
// see e.g. http://www.boost.org/doc/libs/1_66_0/libs/python/doc/html/tutorial/tutorial/techniques.html
//

void def_pymod_common(py::module &mod);

void def_pymod_consts(py::module &mod);
void def_pymod_struct(py::module &mod);

void def_pymod_domain(py::module &mod);

void def_pymod_neighbor(py::module &mod);

void def_pymod_weight(py::module &mod);

void def_pymod_operator(py::module &mod);

//--------------------------------------------------------------------------------


//
// define module
//
PYBIND11_MODULE(pympscore, mod) {

	mod.doc() = "pymps.core";

	////////////////////////////////////////

	def_pymod_common(mod);

	def_pymod_consts(mod);

	def_pymod_struct(mod);

	def_pymod_domain(mod);

	def_pymod_neighbor(mod);

	def_pymod_weight(mod);

	def_pymod_operator(mod);
}

//--------------------------------------------------------------------------------





#pragma warning(pop)
////////////////////////////////////////////////////////////////////////////////



