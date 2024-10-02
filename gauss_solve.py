import ctypes

gauss_library_path = './libgauss.so'

def unpack(A):
    """ Extract L and U parts from A, fill with 0's and 1's """
    n = len(A)
    L = [[A[i][j] for j in range(i)] + [1] + [0 for j in range(i + 1, n)]
         for i in range(n)]

    U = [[0 for j in range(i)] + [A[i][j] for j in range(i, n)]
         for i in range(n)]

    return L, U

def lu_c(A):
    """ Accepts a list of lists A of floats and returns (P, L, U) - the PA=LU decomposition as a tuple. """
    # Load the shared library
    lib = ctypes.CDLL(gauss_library_path)

    n = len(A)
    # Create a 2D array in Python and flatten it
    flat_array_2d = [item for row in A for item in row]

    # Convert to a ctypes array
    c_array_2d = (ctypes.c_double * len(flat_array_2d))(*flat_array_2d)

    # Prepare to retrieve the permutation array
    P = (ctypes.c_int * n)()
    
    # Define the function signature
    lib.plu.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int))

    # Call the C function for LU decomposition
    lib.plu(n, c_array_2d, P)

    # Convert back to a 2D Python list of lists
    modified_array_2d = [
        [c_array_2d[i * n + j] for j in range(n)]
        for i in range(n)
    ]

    # Unpack L and U from the modified array
    L, U = unpack(modified_array_2d)

    # Convert permutation array to a Python list
    P_list = list(P)

    return P_list, L, U

def lu_python(A):
    """ Accepts a list of lists A of floats and returns (P, L, U) - the PA=LU decomposition as a tuple. """
    n = len(A)
    P = list(range(n))  # Initialize the permutation array as an identity permutation
    
    for k in range(n):
        # Find the pivot element
        pivot = max(range(k, n), key=lambda i: abs(A[i][k]))
        if pivot != k:
            A[k], A[pivot] = A[pivot], A[k]  # Swap rows in A
            P[k], P[pivot] = P[pivot], P[k]  # Swap in permutation array

        for i in range(k + 1, n):
            A[i][k] /= A[k][k]  # Calculate the multipliers
            for j in range(k + 1, n):
                A[i][j] -= A[i][k] * A[k][j]  # Eliminate

    return unpack(A), P

def lu(A, use_c=False):
    if use_c:
        return lu_c(A)
    else:
        return lu_python(A)

if __name__ == "__main__":
    def get_A():
        """ Make a test matrix """
        A = [[2.0, 3.0, -1.0],
             [4.0, 1.0, 2.0],
             [-2.0, 7.0, 2.0]]
        return A

    A = get_A()

    P, L, U = lu(A, use_c=False)
    print("Using Python:")
    print("P:", P)
    print("L:", L)
    print("U:", U)

    # Must re-initialize A as it was destroyed
    A = get_A()

    P, L, U = lu(A, use_c=True)
    print("\nUsing C:")
    print("P:", P)
    print("L:", L)
    print("U:", U)
