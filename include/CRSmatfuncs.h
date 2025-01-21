#ifndef CRSMATFUNCS_H
#define CRSMATFUNCS_H
    #include "CRSMat_types.h"
    void csr_matvec(CRSMatrix *A, double *x, double *y);
    void choleskyDecomposition(CRSMatrix *A, CRSMatrix *L);
    void forward_substitution(CRSMatrix *L, double *b, double *y);
    void backward_substitution(CRSMatrix *L, double *y, double *x);
    int isPositiveDefiniteMinor(int k, double *values, int *columns, int *row_ptr);
    int isPositiveDefinite(int n, double *values, int *columns, int *row_ptr);
    int countNonZero(int rows, int cols, double *matrix);
    void convertToCRS(int rows, int cols, double *matrix, double *val, int *col_ind, int *row_ptr);
    void print_CRSMatrix(CRSMatrix *mat) ;
    void transCRSmat (CRSMatrix *A, CRSMatrix *AT);
#endif