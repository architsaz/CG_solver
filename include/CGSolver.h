#ifndef CGSOLVER_H
#define CGSOLVER_H
    #include "types.h"
    void conjugate_gradient(CRSMatrix *A, double *b, double *u);
    void precond_conjugate_gradient(CRSMatrix *A, double *b, double *u);
    void apply_preconditioner(CRSMatrix *A, double *r, double *z);
    void solver_set_config(SolverConfig new_config);
#endif