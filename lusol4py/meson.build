python_sources = [
  '__init__.py',
  'lusol.py',
]

py.install_sources(
  python_sources,
  subdir: 'lusol4py'
)

if enable_fortran
  if blas_dep.found() and lapack_dep.found()
    message('Add fortran wrapper since all dependencies are found')
    subdir('dummy')
    subdir('internal')
  else
    message('Skipp fortran wrapper since not dependencies are found')
  endif
endif