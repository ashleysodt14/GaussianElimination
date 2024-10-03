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
    for (int i = 0; i < n; i++) {
        P[i] = i; // Initialize the permutation matrix
    }

    for (int k = 0; k < n; k++) {
        // Pivot
        int maxIndex = k;
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > fabs(A[maxIndex][k])) {
                maxIndex = i;
            }
        }

        if (maxIndex != k) {
            swap_rows(A, P, k, maxIndex, n);
        }

        // Elimination
        for (int i = k + 1; i < n; i++) {
            A[i][k] /= A[k][k];
            for (int j = k + 1; j < n; j++) {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }
}

void print_matrix(int n, double A[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
}

void print_permutation(int n, int P[n]) {
    for (int i = 0; i < n; i++) {
        printf("%d ", P[i]);
    }
    printf("\n");
}
