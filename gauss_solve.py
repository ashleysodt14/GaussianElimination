#----------------------------------------------------------------
# File:     gauss_solve.py
#----------------------------------------------------------------
#
# Author:   Marek Rychlik (rychlik@arizona.edu)
# Date:     Thu Sep 26 10:38:32 2024
# Copying:  (C) Marek Rychlik, 2020. All rights reserved.
# 
#----------------------------------------------------------------
# A Python wrapper module around the C library libgauss.so

import ctypes

gauss_library_path = './libgauss.so'

def unpack(A):
    """ Extract L and U parts from A, fill with 0's and 1's """
    n = len(A)
    L = [[A[i][j] for j in range(i)] + [1] + [0 for j in range(i+1, n)]
         for i in range(n)]

    U = [[0 for j in range(i)] + [A[i][j] for j in range(i, n)]
         for i in range(n)]

    return L, U

def lu_c(A):
    """ Accepts a list of lists A of floats and
    it returns (L, U) - the LU-decomposition as a tuple.
    """
    # Load the shared library
    lib = ctypes.CDLL(gauss_library_path)

    # Create a 2D array in Python and flatten it
    n = len(A)
    flat_array_2d = [item for row in A for item in row]

    # Convert to a ctypes array
    c_array_2d = (ctypes.c_double * len(flat_array_2d))(*flat_array_2d)

    # Define the function signature
    lib.lu_in_place.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double))

    # Modify the array in C (e.g., add 10 to each element)
    lib.lu_in_place(n, c_array_2d)

    # Convert back to a 2D Python list of lists
    modified_array_2d = [
        [c_array_2d[i * n + j] for j in range(n)]
        for i in range(n)
    ]

    # Extract L and U parts from A, fill with 0's and 1's
    return unpack(modified_array_2d)

def lu_python(A):
    n = len(A)
    for k in range(n):
        for i in range(k,n):
            for j in range(k):
                A[k][i] -= A[k][j] * A[j][i]
        for i in range(k+1, n):
            for j in range(k):
                A[i][k] -= A[i][j] * A[j][k]
            A[i][k] /= A[k][k]

    return unpack(A)


def lu(A, use_c=False):
    if use_c:
        return lu_c(A)
    else:
        return lu_python(A)

def plu_python(A):
    """PA=LU decomposition in pure Python."""
    n = len(A)
    P = list(range(n))
    L = [[0.0] * n for _ in range(n)]
    U = [row[:] for row in A]  # Make a copy of A

    for k in range(n - 1):
        # Find the pivot element
        pivot_row = max(range(k, n), key=lambda i: abs(U[i][k]))
        if U[pivot_row][k] == 0:
            raise ValueError("Matrix is singular and cannot be decomposed.")

        # Swap rows in U and P
        if pivot_row != k:
            U[k], U[pivot_row] = U[pivot_row], U[k]
            P[k], P[pivot_row] = P[pivot_row], P[k]

        # Eliminate below the pivot
        for i in range(k + 1, n):
            L[i][k] = U[i][k] / U[k][k]
            for j in range(k, n):
                U[i][j] -= L[i][k] * U[k][j]

    # Set the diagonal of L to 1
    for i in range(n):
        L[i][i] = 1.0

    return P, L, U

def plu_c(A):
    """PA=LU decomposition using the C implementation."""
    n = len(A)
    
    # Load the shared C library (adjust the path as necessary)
    lib = ctypes.CDLL(gauss_library_path)

    # Convert Python lists to C-compatible data structures
    A_c = (ctypes.POINTER(ctypes.c_double) * n)()
    for i in range(n):
        A_c[i] = (ctypes.c_double * n)(*A[i])

    P_c = (ctypes.c_int * n)()

    # Call the C function
    lib.plu(ctypes.c_int(n), A_c, P_c)

    # Extract the results back to Python lists
    P = [P_c[i] for i in range(n)]
    L = [[0.0] * n for _ in range(n)]
    U = [[A_c[i][j] for j in range(n)] for i in range(n)]
    
    # Construct L (since C stores L's multipliers directly in A)
    for i in range(n):
        L[i][i] = 1.0
        for j in range(i):
            L[i][j] = A_c[i][j]

    return P, L, U

def plu(A, use_c=False):
    """
    PA=LU decomposition function.
    
    Parameters:
    A : list of lists
        The matrix to decompose.
    use_c : bool
        If True, use the C implementation, otherwise use Python implementation.
    
    Returns:
    P, L, U : lists
        Permutation list, lower triangular matrix L, upper triangular matrix U.
    """
    if use_c:
        return plu_c(A)
    else:
        return plu_python(A)


if __name__ == "__main__":

    def get_A():
        """ Make a test matrix """
        A = [[2.0, 3.0, -1.0],
             [4.0, 1.0, 2.0],
             [-2.0, 7.0, 2.0]]
        return A

    A = get_A()

    L, U = lu(A, use_c = False)
    print(L)
    print(U)

    # Must re-initialize A as it was destroyed
    A = get_A()

    L, U = lu(A, use_c=True)
    print(L)
    print(U)
