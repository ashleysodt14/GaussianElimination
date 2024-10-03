#include <stdio.h>
#include <math.h>

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

// Helper function to print a matrix
void print_matrix(int n, double *A) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", A[i * n + j]);
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
