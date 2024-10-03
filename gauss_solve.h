#ifndef GAUSS_SOLVE_H
#define GAUSS_SOLVE_H

// Function to perform PLU decomposition in-place
void plu(int n, double *A, int P[n]);

// Helper function to swap rows in a matrix
void swap_rows(double *A, int n, int P[], int row1, int row2);

// Helper function to print the matrix (for debugging)
void print_matrix(int n, double *A);

// Helper function to print the permutation matrix (for debugging)
void print_permutation(int n, int P[n]);

#endif
