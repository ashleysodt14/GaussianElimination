#ifndef GAUSS_SOLVE_H
#define GAUSS_SOLVE_H

// Function to perform PLU decomposition in-place
void plu(int n, double A[n][n], int P[n]);

// Helper function to swap rows in a matrix
void swap_rows(double A[][n], int P[], int row1, int row2, int n);

// Helper function to print the matrix (for debugging)
void print_matrix(int n, double A[n][n]);

// Helper function to print the permutation matrix (for debugging)
void print_permutation(int n, int P[n]);

#endif
