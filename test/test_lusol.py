from lusol4py import LUSOL
from scipy.sparse import rand

if __name__ == "__main__":
    m = 5
    n = 5
    density = 0.5
    A = rand(m, n, density, format="coo")
    print(f"A:\n{A}")
    print(f"A.toarray():\n{A.toarray()}")
    # print(f"lu1fac(): {lu1fac(A)}")

    lu = LUSOL(A)
    info = lu.factorize()
    print(f"info: {info}")
