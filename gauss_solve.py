import numpy as np
from ctypes import CDLL, POINTER, c_double, c_int

# Load the C library (ensure you compiled the C code to a shared object)
lib = CDLL('./libgauss.so')

# Define the LU decomposition function
def lu(A, use_c=True):
    n = len(A)

    if use_c:
        A_c = (c_double * (n * n))()

        for i in range(n):
            for j in range(n):
                A_c[i * n + j] = A[i][j]

        lib.lu_in_place(c_int(n), A_c)

        L = np.zeros((n, n))
        U = np.zeros((n, n))

        for i in range(n):
            for j in range(n):
                if i > j:
                    L[i][j] = A_c[i * n + j]
                elif i == j:
                    L[i][j] = 1
                    U[i][j] = A_c[i * n + j]
                else:
                    U[i][j] = A_c[i * n + j]

        return L, U
    else:
        A = np.array(A, dtype=float)
        L = np.zeros((n, n))
        U = np.copy(A)

        for k in range(n):
            for i in range(k + 1, n):
                L[i, k] = U[i, k] / U[k, k]
                U[i, k:] -= L[i, k] * U[k, k:]

        np.fill_diagonal(L, 1)
        return L, U

def plu(A, use_c=True):
    n = len(A)

    if use_c:
        A_c = (c_double * (n * n))()
        P_c = (c_int * n)()

        for i in range(n):
            for j in range(n):
                A_c[i * n + j] = A[i][j]

        lib.plu(c_int(n), A_c, P_c)

        P = [P_c[i] for i in range(n)]
        L = np.zeros((n, n))
        U = np.zeros((n, n))

        for i in range(n):
            for j in range(n):
                if i > j:
                    L[i][j] = A_c[i * n + j]
                elif i == j:
                    L[i][j] = 1
                    U[i][j] = A_c[i * n + j]
                else:
                    U[i][j] = A_c[i * n + j]

        return P, L, U
    else:
        A = np.array(A, dtype=float)
        P = np.arange(n)
        L = np.zeros((n, n))
        U = np.copy(A)

        for k in range(n):
            max_index = np.argmax(np.abs(U[k:, k])) + k
            if max_index != k:
                U[[k, max_index]] = U[[max_index, k]]
                P[[k, max_index]] = P[[max_index, k]]

            for i in range(k + 1, n):
                L[i, k] = U[i, k] / U[k, k]
                U[i, k:] -= L[i, k] * U[k, k:]

        np.fill_diagonal(L, 1)
        return P.tolist(), L, U
