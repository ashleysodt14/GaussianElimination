#include <stdio.h>
#include <math.h>

void swap_rows(double A[][n], int P[], int row1, int row2, int n) {
    for (int i = 0; i < n; i++) {
        double temp = A[row1][i];
        A[row1][i] = A[row2][i];
        A[row2][i] = temp;
    }
    int tempP = P[row1];
    P[row1] = P[row2];
    P[row2] = tempP;
}

void plu(int n, double A[n][n], int P[n]) {
    // Initialize permutation matrix P to the identity matrix
    for (int i = 0; i < n; i++) {
        P[i] = i;
    }

    for (int k = 0; k < n; k++) {
        // Find the pivot element
        int maxIndex = k;
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > fabs(A[maxIndex][k])) {
                maxIndex = i;
            }
        }

        // Swap rows in U and P if necessary
        if (maxIndex != k) {
            swap_rows(A, P, k, maxIndex, n);
        }

        // Perform the elimination process
        for (int i = k + 1; i < n; i++) {
            // Compute the multiplier and store it in the lower triangular part of A
            A[i][k] /= A[k][k];

            // Update U (upper triangular part of A)
            for (int j = k + 1; j < n; j++) {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }

    // At this point, A contains both L (in the lower part) and U (in the upper part)
}

// Helper function to print a matrix
void print_matrix(int n, double A[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
}

void print_permutation(int n, int P[]) {
    for (int i = 0; i < n; i++) {
        printf("%d ", P[i]);
    }
    printf("\n");
}
