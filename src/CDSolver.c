#include <stdlib.h>
#include <stdio.h>
#include "CRSMat_types.h"
#include "CRSmatfuncs.h"

// Function to solve Ax = b using Cholesky decomposition
void solve_cholesky(CRSMatrix *A, double *b, double *x) {
    CRSMatrix L;
    choleskyDecomposition(A, &L);
    #ifdef DEBUG
        printf("L:\n");
        print_CRSMatrix(&L);
    #endif

    double *y = (double *)malloc((size_t)A->n * sizeof(double));
    forward_substitution(&L, b, y);

    #ifdef DEBUG
    printf("Solution y:\n");
    for (int i = 0; i < A->n; i++){
        printf("%f ", y[i]);
        printf("\n");
    }
    #endif

    CRSMatrix LT;
    transCRSmat (&L,&LT);

    #ifdef DEBUG
    printf("LT:\n");
    print_CRSMatrix(&LT);
    #endif

    backward_substitution(&LT, y, x);

    free(y);
    free(L.row_ptr);
    free(L.col_index);
    free(L.values);
    free(LT.row_ptr);
    free(LT.col_index);
    free(LT.values);
}

