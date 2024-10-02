#include "gauss_solve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h> // Include math.h for fabs function

void plu(int n, double A[n][n], int P[n]) {
    double L[n][n];  // Temporary lower triangular matrix
    double U[n][n];  // Temporary upper triangular matrix

    // Initialize P with the identity permutation
    for (int i = 0; i < n; i++) {
        P[i] = i;
    }

    // LU decomposition with partial pivoting
    for (int k = 0; k < n; k++) {
        // Find the pivot
        double max = fabs(A[k][k]);
        int maxIndex = k;
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > max) {
                max = fabs(A[i][k]);
                maxIndex = i;
            }
        }

        // Swap rows in A and update the permutation
        if (maxIndex != k) {
            for (int j = 0; j < n; j++) {
                double temp = A[k][j];
                A[k][j] = A[maxIndex][j];
                A[maxIndex][j] = temp;
            }

            // Swap the permutation
            int temp = P[k];
            P[k] = P[maxIndex];
            P[maxIndex] = temp;
        }

        // Decompose into L and U
        for (int i = k; i < n; i++) {
            U[k][i] = A[k][i];  // Upper triangular part
        }
        for (int i = k + 1; i < n; i++) {
            L[i][k] = A[i][k] / U[k][k];  // Lower triangular part
            for (int j = k; j < n; j++) {
                A[i][j] -= L[i][k] * U[k][j];
            }
        }
    }

    // Copy U and L to the appropriate positions in A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i < j) {
                A[i][j] = U[i][j];  // Upper triangular part
            } else if (i == j) {
                A[i][j] = 1.0;  // Diagonal of L is 1
            } else {
                A[i][j] = L[i][j];  // Lower triangular part
            }
        }
    }
}

void gauss_solve_in_place(const int n, double A[n][n], double b[n]) {
    for (int k = 0; k < n; ++k) {
        for (int i = k + 1; i < n; ++i) {
            /* Store the multiplier into A[i][k] as it would become 0 and be
            useless */
            A[i][k] /= A[k][k];
            for (int j = k + 1; j < n; ++j) {
                A[i][j] -= A[i][k] * A[k][j];
            }
            b[i] -= A[i][k] * b[k];
        }
    } /* End of Gaussian elimination, start back-substitution. */
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            b[i] -= A[i][j] * b[j];
        }
        b[i] /= A[i][i];
    } /* End of back-substitution. */
}

void lu_in_place(const int n, double A[n][n]) {
    for (int k = 0; k < n; ++k) {
        for (int i = k; i < n; ++i) {
            for (int j = 0; j < k; ++j) {
                /* U[k][i] -= L[k][j] * U[j][i] */
                A[k][i] -= A[k][j] * A[j][i];
            }
        }
        for (int i = k + 1; i < n; ++i) {
            for (int j = 0; j < k; ++j) {
                /* L[i][k] -= A[i][k] * U[j][k] */
                A[i][k] -= A[i][j] * A[j][k];
            }
            /* L[k][k] /= U[k][k] */
            A[i][k] /= A[k][k];
        }
    }
}

void lu_in_place_reconstruct(int n, double A[n][n]) {
    for (int k = n - 1; k >= 0; --k) {
        for (int i = k + 1; i < n; ++i) {
            A[i][k] *= A[k][k];
            for (int j = 0; j < k; ++j) {
                A[i][k] += A[i][j] * A[j][k];
            }
        }
        for (int i = k; i < n; ++i) {
            for (int j = 0; j < k; ++j) {
                A[k][i] += A[k][j] * A[j][i];
            }
        }
    }
}
