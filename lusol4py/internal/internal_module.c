// #include <Python.h>
// #define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
// #include "numpy/arrayobject.h"

// #ifdef HAVE_BLAS_ILP64
// #define F_INT npy_int64
// #define F_INT_NPY NPY_INT64
// #else
// #define F_INT int
// #define F_INT_NPY NPY_INT
// #endif

// #if defined(UPPERCASE_FORTRAN)
//     #if defined(NO_APPEND_FORTRAN)
//         /* nothing to do here */
//     #else
//         #define LU1FAC LU1FAC_
//     #endif
// #else
//     #if defined(NO_APPEND_FORTRAN)
//         #define LU1FAC lu1fac
//     #else
//         #define LU1FAC lu1fac_
//     #endif
// #endif

// typedef void lu1fac_t(F_INT *m, F_INT *n, F_INT *nelem, F_INT *lena, 
//                       F_INT *luparm, double *paramlu, double *a, F_INT *indc,
//                       F_INT *indr, F_INT *p, F_INT *q, F_INT *lenc, F_INT *lenr,
//                       F_INT *locc, F_INT *locr, F_INT *iploc, F_INT *iqloc, 
//                       F_INT *ipinv, F_INT *iqinv, double *w, F_INT *inform);

// extern lu1fac_t LU1FAC;

// static PyObject* py_lu1fac(PyObject* self, PyObject* args) {
//     PyArrayObject *luparm, *parmlu, *a, *indc, *indr;
//     PyArrayObject *p, *q, *lenc, *lenr, *locc, *locr;
//     PyArrayObject *iploc, *iqloc, *ipinv, *iqinv, *w;

//     int m, n, nelem, lena;

//     // parse Python args
//     if (!PyArg_ParseTuple(args, "iiiiOOOOOOOOOOOOOOOO",
//                           &m, &n, &nelem, &lena,
//                           &luparm, &parmlu,
//                           &a, &indc, &indr,
//                           &p, &q,
//                           &lenc, &lenr, &locc, &locr,
//                           &iploc, &iqloc, &ipinv, &iqinv,
//                           &w)) {
//         return NULL;
//     }

//     // ensure arrays are C-contiguous and correctly typed
//     // #define CHECK_ARRAY(arr, typecode) \
//     //     if (!PyArray_Check(arr) || PyArray_TYPE(arr) != typecode || !PyArray_ISCONTIGUOUS(arr)) { \
//     //         PyErr_Format(PyExc_TypeError, "Expected contiguous array of type %s", PyArray_DescrFromType(typecode)->type); \
//     //         return NULL; \
//     //     }
//     #define CHECK_ARRAY(obj, typecode)                                      \
//     do {                                                                \
//         if (obj == NULL || !PyArray_Check(obj)) {                       \
//             PyErr_SetString(PyExc_TypeError, "Expected numpy array");  \
//             return NULL;                                                \
//         }                                                               \
//         if (PyArray_TYPE(obj) != typecode || !PyArray_ISCONTIGUOUS(obj)) { \
//             PyErr_Format(PyExc_TypeError,                               \
//                          "Expected contiguous array of type %s",       \
//                          PyArray_DescrFromType(typecode)->type);       \
//             return NULL;                                                \
//         }                                                               \
//     } while (0)

//     // CHECK_ARRAY(luparm, NPY_INT64);
//     CHECK_ARRAY(parmlu, NPY_FLOAT64);
//     // CHECK_ARRAY(a, NPY_FLOAT64);
//     // CHECK_ARRAY(indc, NPY_INT64);
//     // CHECK_ARRAY(indr, NPY_INT64);
//     // CHECK_ARRAY(p, NPY_INT64);
//     // CHECK_ARRAY(q, NPY_INT64);
//     // CHECK_ARRAY(lenc, NPY_INT64);
//     // CHECK_ARRAY(lenr, NPY_INT64);
//     // CHECK_ARRAY(locc, NPY_INT64);
//     // CHECK_ARRAY(locr, NPY_INT64);
//     // CHECK_ARRAY(iploc, NPY_INT64);
//     // CHECK_ARRAY(iqloc, NPY_INT64);
//     // CHECK_ARRAY(ipinv, NPY_INT64);
//     // CHECK_ARRAY(iqinv, NPY_INT64);
//     // CHECK_ARRAY(w, NPY_FLOAT64);

//     // call Fortran function
//     int inform = 0;
//     // LU1FAC(&m, &n, &nelem, &lena,
//     //        (int*)PyArray_DATA(luparm),
//     //        (double*)PyArray_DATA(parmlu),
//     //        (double*)PyArray_DATA(a),
//     //        (int*)PyArray_DATA(indc),
//     //        (int*)PyArray_DATA(indr),
//     //        (int*)PyArray_DATA(p),
//     //        (int*)PyArray_DATA(q),
//     //        (int*)PyArray_DATA(lenc),
//     //        (int*)PyArray_DATA(lenr),
//     //        (int*)PyArray_DATA(locc),
//     //        (int*)PyArray_DATA(locr),
//     //        (int*)PyArray_DATA(iploc),
//     //        (int*)PyArray_DATA(iqloc),
//     //        (int*)PyArray_DATA(ipinv),
//     //        (int*)PyArray_DATA(iqinv),
//     //        (double*)PyArray_DATA(w),
//     //        &inform);

//     return PyLong_FromLong(inform);
// }

// static PyMethodDef methods[] = {
//     {"lu1fac", (PyCFunction)py_lu1fac, METH_VARARGS, "Sparse LU factorization"},
//     {NULL, NULL, 0, NULL},
// };

// static struct PyModuleDef module = {
//     PyModuleDef_HEAD_INIT,
//     "lusol4py",
//     NULL,
//     -1,
//     methods,
// };

// PyMODINIT_FUNC PyInit_internal(void)
// {
//     return PyModule_Create(&module);
// }


#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include "numpy/arrayobject.h"

#define F_INT npy_int64
#define F_INT_NPY NPY_INT64

// #if defined(UPPERCASE_FORTRAN)
//     #if defined(NO_APPEND_FORTRAN)
//         /* nothing to do here */
//     #else
//         #define LU1FAC LU1FAC_
//     #endif
// #else
//     #if defined(NO_APPEND_FORTRAN)
//         #define LU1FAC lu1fac
//     #else
//         #define LU1FAC lu1fac_
//     #endif
// #endif
#if defined(UPPERCASE_FORTRAN)
    #if defined(NO_APPEND_FORTRAN)
        #define LU1FAC LU1FAC_C
    #else
        #define LU1FAC LU1FAC_C_
    #endif
#else
    #if defined(NO_APPEND_FORTRAN)
        #define LU1FAC lu1fac_c
    #else
        #define LU1FAC lu1fac_c_
    #endif
#endif

typedef void lu1fac_t(F_INT *m, F_INT *n, F_INT *nelem, F_INT *lena, 
                      F_INT *luparm, double *parmlu, double *a, F_INT *indc,
                      F_INT *indr, F_INT *p, F_INT *q, F_INT *lenc, F_INT *lenr,
                      F_INT *locc, F_INT *locr, F_INT *iploc, F_INT *iqloc, 
                      F_INT *ipinv, F_INT *iqinv, double *w, F_INT *inform);

extern lu1fac_t LU1FAC;

#define CHECK_ARRAY(obj, typecode)                                      \
    do {                                                                \
        if (obj == NULL || !PyArray_Check(obj)) {                       \
            PyErr_SetString(PyExc_TypeError, "Expected numpy array");  \
            return NULL;                                                \
        }                                                               \
        if (PyArray_TYPE((PyArrayObject *)obj) != typecode || !PyArray_ISCONTIGUOUS((PyArrayObject *)obj)) { \
            PyErr_Format(PyExc_TypeError,                               \
                         "Expected contiguous array of type %s",       \
                         PyArray_DescrFromType(typecode)->type);       \
            return NULL;                                                \
        }                                                               \
    } while (0)

static PyObject* py_lu1fac(PyObject* self, PyObject* args) {
    PyArrayObject *luparm, *parmlu, *a, *indc, *indr;
    PyArrayObject *p, *q, *lenc, *lenr, *locc, *locr;
    PyArrayObject *iploc, *iqloc, *ipinv, *iqinv, *w;

    F_INT m, n, nelem, lena;

    if (!PyArg_ParseTuple(args, "llllOOOOOOOOOOOOOOOO",
                          &m, &n, &nelem, &lena,
                          &luparm, &parmlu,
                          &a, &indc, &indr,
                          &p, &q,
                          &lenc, &lenr, &locc, &locr,
                          &iploc, &iqloc, &ipinv, &iqinv,
                          &w)) {
        return NULL;
    }

    // check types and contiguity
    CHECK_ARRAY(luparm, F_INT_NPY);
    CHECK_ARRAY(parmlu, NPY_FLOAT64);
    CHECK_ARRAY(a, NPY_FLOAT64);
    CHECK_ARRAY(indc, F_INT_NPY);
    CHECK_ARRAY(indr, F_INT_NPY);
    CHECK_ARRAY(p, F_INT_NPY);
    CHECK_ARRAY(q, F_INT_NPY);
    CHECK_ARRAY(lenc, F_INT_NPY);
    CHECK_ARRAY(lenr, F_INT_NPY);
    CHECK_ARRAY(locc, F_INT_NPY);
    CHECK_ARRAY(locr, F_INT_NPY);
    CHECK_ARRAY(iploc, F_INT_NPY);
    CHECK_ARRAY(iqloc, F_INT_NPY);
    CHECK_ARRAY(ipinv, F_INT_NPY);
    CHECK_ARRAY(iqinv, F_INT_NPY);
    CHECK_ARRAY(w, NPY_FLOAT64);

    F_INT inform = 0;

    LU1FAC(&m, &n, &nelem, &lena,
           (F_INT*)PyArray_DATA(luparm),
           (double*)PyArray_DATA(parmlu),
           (double*)PyArray_DATA(a),
           (F_INT*)PyArray_DATA(indc),
           (F_INT*)PyArray_DATA(indr),
           (F_INT*)PyArray_DATA(p),
           (F_INT*)PyArray_DATA(q),
           (F_INT*)PyArray_DATA(lenc),
           (F_INT*)PyArray_DATA(lenr),
           (F_INT*)PyArray_DATA(locc),
           (F_INT*)PyArray_DATA(locr),
           (F_INT*)PyArray_DATA(iploc),
           (F_INT*)PyArray_DATA(iqloc),
           (F_INT*)PyArray_DATA(ipinv),
           (F_INT*)PyArray_DATA(iqinv),
           (double*)PyArray_DATA(w),
           &inform);

    return PyLong_FromLong((long)inform);
}

static PyMethodDef methods[] = {
    {"lu1fac", (PyCFunction)py_lu1fac, METH_VARARGS, "Sparse LU factorization"},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "lusol4py",  // module name
    NULL,
    -1,
    methods,
};

PyMODINIT_FUNC PyInit_internal(void)
{
    import_array();  // important for NumPy
    return PyModule_Create(&module);
}
