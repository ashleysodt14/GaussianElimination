import numpy as np

def plu(A):
    n = len(A)
    P = np.arange(n)  # Initialize permutation matrix P as an array of indices
    
    # Perform PLU decomposition
    for k in range(n):
        # Find the pivot element
        max_index = np.argmax(np.abs(A[k:n, k])) + k
        
        # Swap rows in U (A) and P
        if max_index != k:
            A[[k, max_index]] = A[[max_index, k]]
            P[[k, max_index]] = P[[max_index, k]]
        
        # Perform the elimination process
        for i in range(k + 1, n):
            A[i, k] /= A[k, k]  # Store the multiplier in the lower triangular part (L)
            A[i, k+1:] -= A[i, k] * A[k, k+1:]  # Update U (upper triangular part)
    
    # Return P, L (lower part of A), and U (upper part of A)
    L = np.tril(A, -1) + np.eye(n)
    U = np.triu(A)
    
    return P, L, U

# Example usage:
A = np.array([[2.0, 3.0, -1.0],
              [4.0, 1.0, 2.0],
              [-2.0, 7.0, 2.0]])

P, L, U = plu(A)

print("P:", P)
print("L:", L)
print("U:", U)
