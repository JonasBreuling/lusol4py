from warnings import warn
import numpy as np
from scipy.sparse import coo_matrix, coo_array, csr_array
import lusol4py


def _permutation_matrix(p):
    p = np.asarray(p) - 1  # LUSOL uses 1-based indexing
    n = len(p)
    return csr_array((np.ones(n), (np.arange(n), p)), shape=(n, n))

class LUSOL:
    def __init__(self, A, lena=None):
        if not isinstance(A, (coo_matrix, coo_array)):
            warn(
                """Matrix 'A' is not a 'coo_matrix' nor a 'coo_array', we do 
                a possibly conversation here!"""
            )
            A = coo_array(A)

        # dimension, number of nonzeros and estimate of required space
        self.m, self.n = A.shape
        self.nelem = A.nnz
        if lena is None:
            # conservative estimate of required space
            self.lena = max(5 * A.nnz, 10 * (self.m + self.n), A.nnz + 500)
        else:
            self.lena = lena

        # data and index type, see internal/src/lusol_precision.f90
        # dtype = A.dtype
        # itype = A.row.dtype
        dtype = np.float64
        itype = np.int64

        # data, row and colum indices of COO format (fotran counts from 1)
        self.a = np.array(A.data, dtype=dtype)
        self.indc = np.array(A.row, dtype=itype) + 1
        self.indr = np.array(A.col, dtype=itype) + 1

        self.luparm = np.zeros(30, dtype=itype)
        self.luparm[1] = 50 # print level
        self.luparm[2] = 5 # maxcol
        self.luparm[5] = 2 # threshold complete pivoting
        self.luparm[7] = 1 # do not discard factors L and U

        self.parmlu = np.zeros(30, dtype=dtype)
        self.parmlu[0] = 10.0 # max Lij allowed during Factor
        self.parmlu[1] = 10.0 # max Lij allowed during Updates
        self.parmlu[2] = 3.0e-13 # absolute tolerance for treating reals as zero
        self.parmlu[3] = 3.7e-11 # absolute tolerance for flagging small diagonals of U
        self.parmlu[4] = 3.7e-11 # relative tolerance for flagging small diagonals of U
        self.parmlu[5] = 3.0 # Factor limiting waste space in  U
        self.parmlu[6] = 0.3 # density for Markowitz pivot strategy
        self.parmlu[7] = 0.5 # other density for Markowitz pivot strategy

        self.p = np.zeros(self.m, dtype=itype)
        self.q = np.zeros(self.n, dtype=itype)

        self.lenc = np.zeros(self.n, dtype=itype)
        self.lenr = np.zeros(self.m, dtype=itype)
        self.locc = np.zeros(self.n, dtype=itype)
        self.locr = np.zeros(self.m, dtype=itype)
        self.iploc = np.zeros(self.m, dtype=itype)
        self.iqloc = np.zeros(self.n, dtype=itype)
        self.ipinv = np.zeros(self.m, dtype=itype)
        self.iqinv = np.zeros(self.n, dtype=itype)
        self.w = np.zeros(self.n, dtype=dtype)

    def factorize(self):
        info = lusol4py.lu1fac(
            self.m,
            self.n,
            self.nelem,
            self.lena,
            self.luparm,
            self.parmlu,
            self.a,
            self.indc,
            self.indr,
            self.p,
            self.q,
            self.lenc,
            self.lenr,
            self.locc,
            self.locr,
            self.iploc,
            self.iqloc,
            self.ipinv,
            self.iqinv,
            self.w,
        )

        message = {
            0: "LU factors were obtained successfully",
            1: "U appears to be singular, as judged by lu6chk.",
            3: "Some index pair indc(l), indr(l) lies outside the matrix dimensions 1:m , 1:n.",
            4: "Some index pair indc(l), indr(l) duplicates another such pair.",
            7: "The arrays a, indc, indr were not large enough. Their length 'lena' should be increase to at least the value 'minlen' given in luparm(13).",
            8: "There was some other fatal error. (Shouldn't happen!).",
            9: "No diagonal pivot could be found with TSP or TDP. The matrix must not be sufficiently definite or quasi-definite.",
            10: "There was some other fatal error.",
        }

        if info > 0:
            warn(message[info])

        return info
