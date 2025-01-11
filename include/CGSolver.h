#ifndef CGSOLVER_H
#define CGSOLVER_H
    #include "types.h"
    void conjugate_gradient(CSRMatrix *A, double *b, double *u);
    void precond_conjugate_gradient(CSRMatrix *A, double *b, double *u);
    void csr_matvec(CSRMatrix *A, double *x, double *y);
    void apply_preconditioner(CSRMatrix *A, double *r, double *z);
    void convertToCRS(int rows, int cols, double *matrix, double *val, int *col_ind, int *row_ptr);
    int countNonZero(int rows, int cols, double *matrix);
    void solver_set_config(SolverConfig new_config);
#endif