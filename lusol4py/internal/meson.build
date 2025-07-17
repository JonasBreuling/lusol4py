internal_c = [
    'internal_module.c',
]

internal_f = [
    'src/lusol_precision.f90',
    'src/lusol.f90',
    'interface.f90',
]

py.extension_module(
    'internal',
    sources: internal_c + internal_f,
    subdir: 'lusol4py',
    install: true,
    dependencies: [blas_lapack_dep],
)
