#include <stdio.h>
#include "gauss_solve.h"

int main() {
    int n = 3;
    double A[3][3] = {
        {2, -1, 1},
        {3, 3, 9},
        {3, 3, 5}
    };
    int P[3];

    printf("Original Matrix A:\n");
    print_matrix(n, A);

    plu(n, A, P);

    printf("\nMatrix A after PLU decomposition (L and U in place):\n");
    print_matrix(n, A);

    printf("\nPermutation Matrix P:\n");
    print_permutation(n, P);

    return 0;
}
