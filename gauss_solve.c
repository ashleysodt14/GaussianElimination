#include "gauss_solve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // Include math.h for fabs function
#include <Python.h>  // Include Python API

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
            if
