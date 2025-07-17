# see https://numpy.org/doc/stable/f2py/buildtools/meson.html
project(
  'lusol4py', 
  'c', 'fortran',
  default_options: ['warning_level=2']
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

enable_fortran = true

if enable_fortran
  message('enable_fortran=true')
  blas_dep = dependency('blas', required: false)
  lapack_dep = dependency('lapack', required: false)

  # if fortran_compiler.found() and blas_dep.found() and lapack_dep.found()
  if blas_dep.found() and lapack_dep.found()
    message('blas_dep.found() and lapack_dep.found()')
    _args_blas_lapack = []
    blas_lapack_dep = declare_dependency(
      dependencies: [blas_dep, lapack_dep],
      compile_args: _args_blas_lapack
    )
  endif
endif

subdir('lusol4py')