import numpy as np
from ctypes import CDLL, POINTER, c_int, c_double

# Load the C library (make sure this path points to the correct shared library file)
lib = CDLL('./plu_decomposition.so')

def plu(A, use_c=False):
    n = len(A)
    
    if use_c:
        # Call C implementation
        A_c = (c_double * (n * n))()  # Create a 1D array for C
        P_c = (c_int * n)()  # Permutation array
        
        # Populate the C array
        for i in range(n):
            for j in range(n):
                A_c[i * n + j] = A[i][j]
        
        # Call the C function (make sure C function signature matches)
        lib.plu(c_int(n), A_c, P_c)
        
        # Extract the results from the C array
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
        # Python-based PA=LU decomposition
        A = np.array(A, dtype=float)
        P = np.arange(n)
        L = np.zeros((n, n))
        U = np.copy(A)

        for k in range(n):
            # Pivot
            max_index = np.argmax(np.abs(U[k:, k])) + k
            if max_index != k:
                U[[k, max_index]] = U[[max_index, k]]
                P[[k, max_index]] = P[[max_index, k]]

            for i in range(k + 1, n):
                L[i, k] = U[i, k] / U[k, k]
                U[i, k:] -= L[i, k] * U[k, k:]

        np.fill_diagonal(L, 1)
        return P.tolist(), L, U

# Example usage
if __name__ == "__main__":
    A = [[2.0, 3.0, -1.0], [4.0, 1.0, 2.0], [-2.0, 7.0, 2.0]]
    use_c = False
    P, L, U = plu(A, use_c=use_c)
    print("Python-based PLU Decomposition:")
    print("P:", P)
    print("L:", L)
    print("U:", U)
    
    use_c = True
    P, L, U = plu(A, use_c=use_c)
    print("C-based PLU Decomposition:")
    print("P:", P)
    print("L:", L)
    print("U:", U)
