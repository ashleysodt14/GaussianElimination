/*----------------------------------------------------------------
* File:     gauss_solve.c
*----------------------------------------------------------------
*
* Author:   Marek Rychlik (rychlik@arizona.edu)
* Date:     Sun Sep 22 15:40:29 2024
* Copying:  (C) Marek Rychlik, 2020. All rights reserved.
*
*----------------------------------------------------------------*/
#include "gauss_solve.h"
#include <stdio.h>
#include <math.h>

void gauss_solve_in_place(const int n, double A[n][n], double b[n])
{
  for(int k = 0; k < n; ++k) {
    for(int i = k+1; i < n; ++i) {
      /* Store the multiplier into A[i][k] as it would become 0 and be
	 useless */
      A[i][k] /= A[k][k];
      for( int j = k+1; j < n; ++j) {
	A[i][j] -= A[i][k] * A[k][j];
      }
      b[i] -= A[i][k] * b[k];
    }
  } /* End of Gaussian elimination, start back-substitution. */
  for(int i = n-1; i >= 0; --i) {
    for(int j = i+1; j<n; ++j) {
      b[i] -= A[i][j] * b[j];
    }
    b[i] /= A[i][i];
  } /* End of back-substitution. */
}

void lu_in_place(const int n, double A[n][n])
{
  for(int k = 0; k < n; ++k) {
    for(int i = k; i < n; ++i) {
      for(int j=0; j<k; ++j) {
	/* U[k][i] -= L[k][j] * U[j][i] */
	A[k][i] -=  A[k][j] * A[j][i]; 
      }
    }
    for(int i = k+1; i<n; ++i) {
      for(int j=0; j<k; ++j) {
	/* L[i][k] -= A[i][k] * U[j][k] */
	A[i][k] -= A[i][j]*A[j][k]; 
      }
      /* L[k][k] /= U[k][k] */
      A[i][k] /= A[k][k];	
    }
  }
}

void lu_in_place_reconstruct(int n, double A[n][n])
{
  for(int k = n-1; k >= 0; --k) {
    for(int i = k+1; i<n; ++i) {
      A[i][k] *= A[k][k];
      for(int j=0; j<k; ++j) {
	A[i][k] += A[i][j]*A[j][k];
      }
    }
    for(int i = k; i < n; ++i) {
      for(int j=0; j<k; ++j) {
	A[k][i] +=  A[k][j] * A[j][i];
      }
    }
  }
}

#include <math.h>

#include <math.h>
#include <stdbool.h>

#include <math.h>
#include <stdio.h>

void swap_rows(double A[], int n, int row1, int row2) {
    for (int i = 0; i < n; i++) {
        double temp = A[row1 * n + i];
        A[row1 * n + i] = A[row2 * n + i];
        A[row2 * n + i] = temp;
    }
}

void plu(int n, double A[n][n], int P[n]) {
    // Step 1: Initialize permutation vector P to identity
    for (int i = 0; i < n; i++) {
        P[i] = i;
    }

    // PLU decomposition
    for (int k = 0; k < n - 1; k++) {
        // Step 2: Find the pivot element
        int pivot_row = k;
        double max_val = fabs(A[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > max_val) {
                max_val = fabs(A[i][k]);
                pivot_row = i;
            }
        }

        // Step 6: Swap rows in A if needed
        if (pivot_row != k) {
            // Swap rows in A
            for (int j = 0; j < n; j++) {
                SWAP(A[k][j], A[pivot_row][j], double);
            }

            // Step 7: Swap rows in P (which tracks row permutations)
            SWAP(P[k], P[pivot_row], int);
        }

        // Step 9: Perform elimination for rows below the pivot row
        for (int i = k + 1; i < n; i++) {
            // Step 10: Compute the multiplier (store in A below the diagonal for L)
            if (A[k][k] != 0) {
                A[i][k] = A[i][k] / A[k][k];  // Store multiplier in A (L part)
            } else {
                A[i][k] = 0; // Avoid division by zero
            }

            // Step 11-13: Update the matrix A (U part)
            for (int j = k + 1; j < n; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];  // Modify A (U part)
            }
        }
    }

    // Note: L is stored in the lower part of A (below the diagonal),
    // U is stored in the upper part of A (including the diagonal).
}
