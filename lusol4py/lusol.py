from warnings import warn
import numpy as np
from scipy.sparse import coo_matrix, coo_array
import lusol4py  # your C-extension module


class LUSOL:
    def __init__(self, A, lena=None):
        if not isinstance(A, (coo_matrix, coo_array)):
            warn(
                """Matrix 'A' is not a 'coo_matrix' nor a 'coo_array', we do 
                a possibly conversation here!"""
            )
            A = coo_array(A)
        self.A = A  # TODO: is this required?

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

        #  luparm input parameters:                                Typical value
        #
        #  luparm( 1) = nout     File number for printed messages.         6
        #
        #  luparm( 2) = lprint   Print level.                              0
        #                   <  0 suppresses output.
        #                   =  0 gives error messages.
        #                  >= 10 gives statistics about the LU factors.
        #                  >= 50 gives debug output from lu1fac
        #                        (the pivot row and column and the
        #                        no. of rows and columns involved at
        #                        each elimination step).
        #
        #  luparm( 3) = maxcol   lu1fac: maximum number of columns         5
        #                        searched allowed in a Markowitz-type
        #                        search for the next pivot element.
        #                        For some of the factorization, the
        #                        number of rows searched is
        #                        maxrow = maxcol - 1.
        #
        #  luparm( 6) = 0    =>  TPP: Threshold Partial   Pivoting.        0
        #             = 1    =>  TRP: Threshold Rook      Pivoting.
        #             = 2    =>  TCP: Threshold Complete  Pivoting.
        #             = 3    =>  TSP: Threshold Symmetric Pivoting.
        #             = 4    =>  TDP: Threshold Diagonal  Pivoting.
        #                             (TDP not yet implemented).
        #                        TRP and TCP are more expensive than TPP but
        #                        more stable and better at revealing rank.
        #                        Take care with setting parmlu(1), especially
        #                        with TCP.
        #                        NOTE: TSP and TDP are for symmetric matrices
        #                        that are either definite or quasi-definite.
        #                        TSP is effectively TRP for symmetric matrices.
        #                        TDP is effectively TCP for symmetric matrices.
        #
        #  luparm( 8) = keepLU   lu1fac: keepLU = 1 means the numerical    1
        #                        factors will be computed if possible.
        #                        keepLU = 0 means L and U will be discarded
        #                        but other information such as the row and
        #                        column permutations will be returned.
        #                        The latter option requires less storage.

        #  luparm output parameters:
        #
        #  luparm(10) = inform   Return code from last call to any LU routine.
        #  luparm(11) = nsing    No. of singularities marked in the
        #                        output array w(*).
        #  luparm(12) = jsing    Column index of last singularity.
        #  luparm(13) = minlen   Minimum recommended value for  lena.
        #  luparm(14) = maxlen   ?
        #  luparm(15) = nupdat   No. of updates performed by the lu8 routines.
        #  luparm(16) = nrank    No. of nonempty rows of U.
        #  luparm(17) = ndens1   No. of columns remaining when the density of
        #                        the matrix being factorized reached dens1.
        #  luparm(18) = ndens2   No. of columns remaining when the density of
        #                        the matrix being factorized reached dens2.
        #  luparm(19) = jumin    The column index associated with DUmin.
        #  luparm(20) = numL0    No. of columns in initial  L.
        #  luparm(21) = lenL0    Size of initial  L  (no. of nonzeros).
        #  luparm(22) = lenU0    Size of initial  U.
        #  luparm(23) = lenL     Size of current  L.
        #  luparm(24) = lenU     Size of current  U.
        #  luparm(25) = lrow     Length of row file.
        #  luparm(26) = ncp      No. of compressions of LU data structures.
        #  luparm(27) = mersum   lu1fac: sum of Markowitz merit counts.
        #  luparm(28) = nUtri    lu1fac: triangular rows in U.
        #  luparm(29) = nLtri    lu1fac: triangular rows in L.
        #  luparm(30) = nslack   lu1fac: no. of unit vectors at start of U. (info only)
        self.luparm = np.zeros(30, dtype=itype)

        #  parmlu( 1) = Ltol1    Max Lij allowed during Factor.
        #                                                  TPP     10.0 or 100.0
        #                                                  TRP      4.0 or  10.0
        #                                                  TCP      5.0 or  10.0
        #                                                  TSP      4.0 or  10.0
        #                        With TRP and TCP (Rook and Complete Pivoting),
        #                        values less than 25.0 may be expensive
        #                        on badly scaled data.  However,
        #                        values less than 10.0 may be needed
        #                        to obtain a reliable rank-revealing
        #                        factorization.
        #  parmlu( 2) = Ltol2    Max Lij allowed during Updates.            10.0
        #                        during updates.
        #  parmlu( 3) = small    Absolute tolerance for       eps**0.8 = 3.0d-13
        #                        treating reals as zero.
        #  parmlu( 4) = Utol1    Absolute tol for flagging    eps**0.67= 3.7d-11
        #                        small diagonals of U.
        #  parmlu( 5) = Utol2    Relative tol for flagging    eps**0.67= 3.7d-11
        #                        small diagonals of U.
        #                        (eps = machine precision)
        #  parmlu( 6) = Uspace   Factor limiting waste space in  U.      3.0
        #                        In lu1fac, the row or column lists
        #                        are compressed if their length
        #                        exceeds Uspace times the length of
        #                        either file after the last compression.
        #  parmlu( 7) = dens1    The density at which the Markowitz      0.3
        #                        pivot strategy should search maxcol
        #                        columns and no rows.
        #                        (Use 0.3 unless you are experimenting
        #                        with the pivot strategy.)
        #  parmlu( 8) = dens2    the density at which the Markowitz      0.5
        #                        strategy should search only 1 column,
        #                        or (if storage is available)
        #                        the density at which all remaining
        #                        rows and columns will be processed
        #                        by a dense LU code.
        #                        For example, if dens2 = 0.1 and lena is
        #                        large enough, a dense LU will be used
        #                        once more than 10 per cent of the
        #                        remaining matrix is nonzero.

        #  parmlu output parameters:
        #
        #  parmlu(10) = Amax     Maximum element in  A.
        #  parmlu(11) = Lmax     Maximum multiplier in current  L.
        #  parmlu(12) = Umax     Maximum element in current  U.
        #  parmlu(13) = DUmax    Maximum diagonal in  U.
        #  parmlu(14) = DUmin    Minimum diagonal in  U.
        #  parmlu(15) = Akmax    Maximum element generated at any stage
        #                        during TCP factorization.
        #  parmlu(16) = growth   TPP: Umax/Amax    TRP, TCP, TSP: Akmax/Amax
        #  parmlu(17) =
        #  parmlu(18) =
        #  parmlu(19) =
        #  parmlu(20) = resid    lu6sol: residual after solve with U or U'.
        #  ...
        #  parmlu(30) =
        self.parmlu = np.zeros(30, dtype=dtype)

        self.p = np.empty(self.m, dtype=itype)
        self.q = np.empty(self.n, dtype=itype)

        self.lenc = np.empty(self.n, dtype=itype)
        self.lenr = np.empty(self.m, dtype=itype)
        self.locc = np.empty(self.n, dtype=itype)
        self.locr = np.empty(self.m, dtype=itype)
        self.iploc = np.empty(self.n, dtype=itype)
        self.iqloc = np.empty(self.m, dtype=itype)
        self.ipinv = np.empty(self.m, dtype=itype)
        self.iqinv = np.empty(self.n, dtype=itype)
        self.w = np.empty(self.n, dtype=dtype)

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

        #  inform = 0 if the LU factors were obtained successfully.
        #         = 1 if U appears to be singular, as judged by lu6chk.
        #         = 3 if some index pair indc(l), indr(l) lies outside
        #             the matrix dimensions 1:m , 1:n.
        #         = 4 if some index pair indc(l), indr(l) duplicates
        #             another such pair.
        #         = 7 if the arrays a, indc, indr were not large enough.
        #             Their length "lena" should be increase to at least
        #             the value "minlen" given in luparm(13).
        #         = 8 if there was some other fatal error.  (Shouldn't happen!)
        #         = 9 if no diagonal pivot could be found with TSP or TDP.
        #             The matrix must not be sufficiently definite
        #             or quasi-definite.
        #         =10 if there was some other fatal error.

        return info
