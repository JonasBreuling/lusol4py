#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include "numpy/arrayobject.h"

// #include "dassl/dassl.h"
// #include "pside/pside.h"
// #include "radau/radau.h"

// PyDoc_STRVAR(doc,
// "Function signature: solve(F, t_span, y0, yp0, rtol=1e-6, atol=1e-3, J=None, t_eval=None)\n"
// "\n"
// "Solve a DAE system F(t, y, y') = 0.\n"
// "\n"
// "Parameters\n"
// "----------\n"
// "F : callable\n"
// "    A Python function that defines the DAE system. The function \n"
// "    must have the signature `F(t, y, yp)`, where `t` is the\n"
// "    current time, `y` is the state vector, and `yp` is the \n"
// "    derivative of the state vector.\n"
// "\n"
// "t_span : array-like\n"
// "    A 2-element list or array defining the time interval `[t_start, t_end]`\n"
// "    over which to integrate the system.\n"
// "\n"
// "y0 : array-like\n"
// "    The initial conditions for the state vector `y` at the start of the integration.\n"
// "\n"
// "yp0 : array-like\n"
// "    The initial conditions for the derivative of the state vector `yp` at the start of the integration.\n"
// "\n"
// "rtol, atol: float (optional)\n"
// "    The used relative and absolute tolerances. Default values: rtol=1e-6, atol=1e-3.\n"
// "\n"
// "t_eval: array-like (optional)\n"
// "      The requested evaluation points. If not given, 500 equidistance points in t_span are chosen.\n"
// "\n"
// "Returns\n"
// "-------\n"
// "result : dict\n"
// "    A dictionary containing the results of the integration. The dictionary has the following keys:\n"
// "    - 'success': Was the integration successful?\n"
// "    - 'order': List of used integration order.\n"
// "    - 't': List of time points.\n"
// "    - 'y': List of stage values corresponding to t.\n"
// "    - 'yp': List of stage derivatives corresponding to t.\n"
// "    - 'nsteps': Total number of steps.\n"
// "    - 'nf': Number of function evaluations.\n"
// "    - 'njac': Number of jacobian evaluations.\n"
// "    - 'nrejerror': Number of error tests failures.\n"
// "    - 'nrejnewton': Number of convergence tests failures.\n"
// "\n"
// "Examples\n"
// "--------\n"
// "def F(t, y, yp):\n"
// "    # Example: A simple harmonic oscillator\n"
// "    return np.array([\n"
// "        yp[0] - y[1],\n"
// "        yp[1] + y[0],\n"
// "    )]\n"
// "result = integrate(f, [0, 10], [1.0, 0.0], [0.0, -1.0])\n"
// "print(result)\n"
// "print(result['t'])\n"
// "print(result['y'])\n"
// "print(result['yp'])"
// );

// static PyMethodDef methods[] = {
//     {"dassl", (PyCFunction)dassl, METH_VARARGS | METH_KEYWORDS, doc},
//     {"pside", (PyCFunction)pside, METH_VARARGS | METH_KEYWORDS, doc},
//     {"radau5", (PyCFunction)radau5, METH_VARARGS | METH_KEYWORDS, doc},
//     {"radau", (PyCFunction)radau, METH_VARARGS | METH_KEYWORDS, doc},
//     {NULL, NULL, 0, NULL},
// };

// static struct PyModuleDef module = {
//     PyModuleDef_HEAD_INIT,
//     "fortran",
//     NULL,
//     -1,
//     methods,
// };

// PyMODINIT_FUNC PyInit_fortran(void)
// {
//     import_array();
//     return PyModule_Create(&module);
// }


// see https://github.com/scipy/scipy/blob/main/scipy/integrate/_odepackmodule.c
// how a fortran function can be wrapped via a C module that is compiled via meson
#include <Python.h>

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
        #define LU1FAC LU1FAC_
    #endif
#else
    #if defined(NO_APPEND_FORTRAN)
        #define BAR  bar
        #define LU1FAC lu1fac
    #else
        #define BAR  bar_
        #define LU1FAC lu1fac_
    #endif
#endif

extern void BAR();

static PyObject* foo(PyObject* self)
{
    BAR();
    return PyUnicode_FromString("foo");
}

//   subroutine lu1fac( m    , n    , nelem, lena , luparm, parmlu,       &
//                      a    , indc , indr , p    , q     ,               &
//                      lenc , lenr , locc , locr ,                       &
//                      iploc, iqloc, ipinv, iqinv, w     , inform )

//     integer(ip),   intent(in)    :: m, n, nelem, lena

//     integer(ip),   intent(inout) :: luparm(30)
//     integer(ip),   intent(inout) :: indc(lena), indr(lena),            &
//                                     p(m)      , q(n)      ,            &
//                                     lenc(n)   , lenr(m)   ,            &
//                                     iploc(n)  , iqloc(m)  ,            &
//                                     ipinv(m)  , iqinv(n)  ,            &
//                                     locc(n)   , locr(m)
//     real(rp),      intent(inout) :: parmlu(30), a(lena), w(n)

//     integer(ip),   intent(out)   :: inform

typedef void lu1fac_t(F_INT *m, F_INT *n, F_INT *nelem, F_INT *lena, 
                      F_INT *luparm, double *paramlu, double *a, F_INT *indc,
                      F_INT *indr, F_INT *p, F_INT *q, F_INT *lenc, F_INT *lenr,
                      F_INT *locc, F_INT *locr, F_INT *iploc, F_INT *iqloc, 
                      F_INT *ipinv, F_INT *iqinv, double *w, F_INT *inform);

extern lu1fac_t LU1FAC;


// TODO: extract the corresponding values from a scipy sparse matrix/array 
// and only use the extracted numpy arrays here, otherwise we get crazy using 
// python code inside C...
static PyObject* lu1fac(PyObject* self)
{
    // TODO: Do some computations/extractions here.
    return PyUnicode_FromString("foo");
}

static PyMethodDef methods[] = {
    {"foo", (PyCFunction)foo, METH_NOARGS, NULL},
    {"lu1fac", (PyCFunction)lu1fac, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "lusol4py",
    NULL,
    -1,
    methods,
};

PyMODINIT_FUNC PyInit_lusol(void)
{
    return PyModule_Create(&module);
}