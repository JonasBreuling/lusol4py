lusol_src = [
    'src/lusol_precision.f90',
    'src/lusol.f90',
]

dummy_src = [
    'dummy_src/dummy_bar.f90',
]

c_src = [
    'lusol_module.c',
]

py.extension_module(
    'lusol',
    sources: c_src + dummy_src + lusol_src,
    subdir: 'lusol4py',
    install: true,
    dependencies: [blas_lapack_dep],
)
