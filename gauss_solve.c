#include <stdio.h>
#include <math.h>

// Helper to swap rows
void swap_rows(double *A, int n, int P[], int row1, int row2) {
    for (int i = 0; i < n; i++) {
        double temp = A[row1 * n + i];
        A[row1 * n + i] = A[row2 * n + i];
        A[row2 * n + i] = temp;
    }
    int tempP = P[row1];
    P[row1] = P[row2];
    P[row2] = tempP;
}

// PLU decomposition in-place
void plu(int n, double *A, int P[]) {
    for (int i = 0; i < n; i++) {
        P[i] = i;  // Initialize P as the identity matrix
    }

    for (int k = 0; k < n; k++) {
        // Find the pivot element
        int maxIndex = k;
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i * n + k]) > fabs(A[maxIndex * n + k])) {
                maxIndex = i;
            }
        }

        if (maxIndex != k) {
            swap_rows(A, n, P, k, maxIndex);
        }

        for (int i = k + 1; i < n; i++) {
            A[i * n + k] /= A[k * n + k];

            for (int j = k + 1; j < n; j++) {
                A[i * n + j] -= A[i * n + k] * A[k * n + j];
            }
        }
    }
}

// Solve a system of linear equations in-place
void gauss_solve_in_place(int n, double *A, double *b, int P[]) {
    // Perform PLU decomposition in place
    plu(n, A, P);

    // Forward substitution for Ly = Pb
    for (int i = 0; i < n; i++) {
        double sum = b[P[i]];
        for (int j = 0; j < i; j++) {
            sum -= A[i * n + j] * b[j];
        }
        b[i] = sum;
    }

    // Backward substitution for Ux = y
    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[i * n + j] * b[j];
        }
        b[i] = sum / A[i * n + i];
    }
}

// Perform LU decomposition in place
void lu_in_place(int n, double *A) {
    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            A[i * n + k] /= A[k * n + k];

            for (int j = k + 1; j < n; j++) {
                A[i * n + j] -= A[i * n + k] * A[k * n + j];
            }
        }
    }
}

