import numpy as np
from ctypes import CDLL, POINTER, c_double, c_int

# Load the C library (ensure you compiled the C code to a shared object)
lib = CDLL('./libgauss.so')

def plu(A, use_c=True):
    n = len(A)

    # Call C implementation
    if use_c:
        # Convert Python arrays to C-compatible format
        A_c = (c_double * (n * n))()
        P_c = (c_int * n)()

        # Flatten and populate the A matrix for C
        for i in range(n):
            for j in range(n):
                A_c[i * n + j] = A[i][j]

        # Call the C function plu
        lib.plu(c_int(n), A_c, P_c)

        # Convert C output back to Python
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

# Example usage:
A = np.array([[2.0, 3.0, -1.0],
              [4.0, 1.0, 2.0],
              [-2.0, 7.0, 2.0]])

P, L, U = plu(A, use_c=True)

print("P:", P)
print("L:", L)
print("U:", U)
