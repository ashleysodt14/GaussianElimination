import ctypes

# Path to the C shared library
lib_path = './libgauss.so'

def extract_lu(A):
    """ Unpacks L and U matrices from matrix A. """
    n = len(A)
    L = [[A[i][j] if j < i else (1 if i == j else 0) for j in range(n)] for i in range(n)]
    U = [[A[i][j] if j >= i else 0 for j in range(n)] for i in range(n)]
    return L, U

def lu_c_decomposition(A):
    """ Uses the C library to perform LU decomposition. """
    n = len(A)
    flat_A = [val for row in A for val in row]
    c_A = (ctypes.c_double * len(flat_A))(*flat_A)
    
    lib = ctypes.CDLL(lib_path)
    lib.lu_in_place.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double))
    lib.lu_in_place(n, c_A)
    
    new_A = [[c_A[i * n + j] for j in range(n)] for i in range(n)]
    return extract_lu(new_A)

def lu_python_decomposition(A):
    """ LU decomposition in pure Python. """
    n = len(A)
    for k in range(n):
        for i in range(k, n):
            for j in range(k):
                A[k][i] -= A[k][j] * A[j][i]
        for i in range(k + 1, n):
            for j in range(k):
                A[i][k] -= A[i][j] * A[j][k]
            A[i][k] /= A[k][k]
    return extract_lu(A)

def lu(A, use_c=False):
    """ Chooses between Python and C implementation of LU decomposition. """
    if use_c:
        return lu_c_decomposition(A)
    return lu_python_decomposition(A)

def plu_python_decomposition(A):
    """ Performs PA=LU decomposition in Python. """
    n = len(A)
    P = list(range(n))
    L = [[0.0] * n for _ in range(n)]
    U = [row[:] for row in A]

    for k in range(n - 1):
        pivot_row = max(range(k, n), key=lambda i: abs(U[i][k]))
        if U[pivot_row][k] == 0:
            raise ValueError("Singular matrix.")

        if pivot_row != k:
            U[k], U[pivot_row] = U[pivot_row], U[k]
            P[k], P[pivot_row] = P[pivot_row], P[k]
            for j in range(k):
                L[k][j], L[pivot_row][j] = L[pivot_row][j], L[k][j]

        for i in range(k + 1, n):
            L[i][k] = U[i][k] / U[k][k]
            for j in range(k, n):
                U[i][j] -= L[i][k] * U[k][j]

    for i in range(n):
        L[i][i] = 1.0

    return P, L, U

def plu_c_decomposition(A):
    """ Uses the C library to perform PA=LU decomposition. """
    n = len(A)
    flat_A = [val for row in A for val in row]
    P = list(range(n))

    c_A = (ctypes.c_double * len(flat_A))(*flat_A)
    c_P = (ctypes.c_int * n)(*P)

    lib = ctypes.CDLL(lib_path)
    lib.plu.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int))
    lib.plu(n, c_A, c_P)

    new_A = [[c_A[i * n + j] for j in range(n)] for i in range(n)]
    permutation_vector = [int(c_P[i]) for i in range(n)]
    return permutation_vector, *extract_lu(new_A)

def plu(A, use_c=False):
    """ Chooses between Python and C implementation of PLU decomposition. """
    if use_c:
        return plu_c_decomposition(A)
    return plu_python_decomposition(A)

if __name__ == "__main__":
    def create_test_matrix():
        return [[2.0, 3.0, -1.0], [4.0, 1.0, 2.0], [-2.0, 7.0, 2.0]]

    A = create_test_matrix()

    L, U = perform_lu(A, use_c=False)
    print("LU (Python):", L, U)

    A = create_test_matrix()
    L, U = perform_lu(A, use_c=True)
    print("LU (C):", L, U)
