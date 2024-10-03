import numpy as np
from ctypes import CDLL, c_int, POINTER, c_double

# Load the C library if needed (modify the path based on your C library's location)
lib = CDLL('./plu_decomposition.so')

def plu(A, use_c=False):
    n = len(A)
    
    if use_c:
        # Convert Python matrix to C-friendly format
        A_c = (c_double * n * n)()
        P_c = (c_int * n)()
        
        # Fill A_c with values from A
        for i in range(n):
            for j in range(n):
                A_c[i][j] = A[i][j]
        
        # Call the C function
        lib.plu(c_int(n), A_c, P_c)
        
        # Convert C output back to Python
        P = [P_c[i] for i in range(n)]
        L = np.zeros((n, n))
        U = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i > j:
                    L[i][j] = A_c[i][j]
                elif i == j:
                    L[i][j] = 1
                    U[i][j] = A_c[i][j]
                else:
                    U[i][j] = A_c[i][j]
        
        return P, L, U
    
    else:
        # Implement the PLU decomposition in Python
        A = np.array(A, dtype=float)
        P = np.arange(n)
        L = np.zeros((n, n))
        U = np.copy(A)

        for k in range(n):
            # Pivoting
            max_index = np.argmax(abs(U[k:, k])) + k
            if max_index != k:
                U[[k, max_index]] = U[[max_index, k]]
                P[[k, max_index]] = P[[max_index, k]]

            # Elimination
            for i in range(k+1, n):
                L[i, k] = U[i, k] / U[k, k]
                U[i, k:] -= L[i, k] * U[k, k:]

        np.fill_diagonal(L, 1)
        return P.tolist(), L, U

# Example usage
A = [[2.0, 3.0, -1.0], [4.0, 1.0, 2.0], [-2.0, 7.0, 2.0]]
use_c = False
P, L, U = plu(A, use_c=use_c)
print("P:", P)
print("L:", L)
print("U:", U)
