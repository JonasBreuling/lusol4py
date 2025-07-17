#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include "numpy/arrayobject.h"

#ifdef HAVE_BLAS_ILP64
#define F_INT npy_int64
#define F_INT_NPY NPY_INT64
#else
#define F_INT int
#define F_INT_NPY NPY_INT
#endif

#if defined(UPPERCASE_FORTRAN)
    #if defined(NO_APPEND_FORTRAN)
        /* nothing to do here */
    #else
        #define BAR  BAR_
    #endif
#else
    #if defined(NO_APPEND_FORTRAN)
        #define BAR  bar
    #else
        #define BAR  bar_
    #endif
#endif

extern void BAR();

static PyObject* foo(PyObject* self)
{
    BAR();
    return PyUnicode_FromString("foo");
}

static PyMethodDef methods[] = {
    {"foo", (PyCFunction)foo, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "lusol4py",
    NULL,
    -1,
    methods,
};

PyMODINIT_FUNC PyInit_dummy(void)
{
    return PyModule_Create(&module);
}