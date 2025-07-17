from lusol4py import LUSOL
from scipy.sparse import rand, coo_array
import numpy as np

def test_3x3():
    # fmt: off
    dense = np.array([
        [0, 5, 22 / 3],
        [4, 2, 1],
        [2, 7, 9],
    ])
    # fmt: on
    print(dense)
    A = coo_array(dense)
    print(A)

    lu = LUSOL(A)
    info = lu.factorize()
    print(f"info: {info}")
    

def test_random():
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

if __name__ == "__main__":
    test_random()
    # test_3x3()
