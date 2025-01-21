#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "CRSMat_types.h"
#include "CGSolver_types.h"
#include "RTgnuplot.h"
#include "CRSmatfuncs.h"

// Default configuration
static SolverConfig solver_config = {1000, 1e-9, false}; // Default: 1000 iterations, residual 1e-9
// Set configuration
void solver_set_config(SolverConfig new_config) {
    solver_config = new_config;
}
// Define a preconditioner application function:
void apply_preconditioner(CRSMatrix *A, double *r, double *z)
{
    // Jacobi Preconditioner
    for (int i = 0; i < A->n; i++) {
        // Find the diagonal element in row i
        double diag = 0.0;
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
            if (A->col_index[j] == i) { // Check if this is the diagonal entry
                diag = A->values[j];
                break;
            }
        }
        if (diag == 0.0) {
            fprintf(stderr,"Error: Zero diagonal element at row %d\n", i);
            exit(EXIT_FAILURE);
        }
        z[i] = r[i] / diag; // Apply diagonal scaling
    }
}
// Conjugate Gradient Solver with preconditioner
void precond_conjugate_gradient(CRSMatrix *A, double *b, double *u)
{
    int n = A->n;
    double *r = malloc((size_t)n * sizeof(double));
    double *p = malloc((size_t)n * sizeof(double));
    double *Ap = malloc((size_t)n * sizeof(double));
    double *z = malloc((size_t)n * sizeof(double));

    // Initialize
    csr_matvec(A, u, r); // r = Au
    for (int i = 0; i < n; i++)
        r[i] = b[i] - r[i];
    apply_preconditioner(A, r, z);
    for (int i = 0; i < n; i++)
        p[i] = z[i];

    double rs_old = 0.0;
    for (int i = 0; i < n; i++)
        rs_old += r[i] * z[i];
    int k;
    for (k = 0; k < solver_config.max_iteration; k++)
    {
        csr_matvec(A, p, Ap);
        double alpha = rs_old;
        double temp = 0.0;
        for (int i = 0; i < n; i++)
            temp += p[i] * Ap[i];
        alpha /= temp;

        for (int i = 0; i < n; i++)
            u[i] += alpha * p[i];
        for (int i = 0; i < n; i++)
            r[i] -= alpha * Ap[i];

        apply_preconditioner(A, r, z);
        double rs_new = 0.0;
        for (int i = 0; i < n; i++)
            rs_new += r[i] * z[i];
        if (sqrt(rs_new) < solver_config.residual_limit)
        {
            printf("* CG solver converged in %d iterations.\n", k + 1);
            break;
        }

        double beta = rs_new / rs_old;
        for (int i = 0; i < n; i++)
            p[i] = z[i] + beta * p[i];
        rs_old = rs_new;
    }
    if (k == solver_config.max_iteration)
        printf("! CG solver reach to the max of iteration (%d).\n", solver_config.max_iteration);

    free(r);
    free(p);
    free(Ap);
    free(z);
}
// Conjugate Gradient Solver
void conjugate_gradient(CRSMatrix *A, double *b, double *u)
{
    int n = A->n;
    double *r = malloc((size_t)n * sizeof(double));
    double *p = malloc((size_t)n * sizeof(double));
    double *Ap = malloc((size_t)n * sizeof(double));
    // Initialize function pointers
    FILE *gnuplot = NULL;
    void (*update_plot)(FILE *, int, double) = noop_update;
    void (*finalize_plot)(FILE *) = noop_finalize;
    if(solver_config.showplot){
        // Initialize Gnuplot
        gnuplot = initialize_gnuplot();
        update_plot = update_gnuplot;
        finalize_plot = finalize_gnuplot;
    }
    // Initialize
    csr_matvec(A, u, r); // r = Au
    for (int i = 0; i < n; i++)
        r[i] = b[i] - r[i];
    for (int i = 0; i < n; i++)
        p[i] = r[i];

    double rs_old = 0.0;
    for (int i = 0; i < n; i++)
        rs_old += r[i] * r[i];
    int k;
    for (k = 0; k < solver_config.max_iteration; k++)
    {
        csr_matvec(A, p, Ap);
        double alpha = rs_old;
        double temp = 0.0;
        for (int i = 0; i < n; i++)
            temp += p[i] * Ap[i];
        alpha /= temp;

        for (int i = 0; i < n; i++)
            u[i] += alpha * p[i];
        for (int i = 0; i < n; i++)
            r[i] -= alpha * Ap[i];

        double rs_new = 0.0;
        for (int i = 0; i < n; i++)
            rs_new += r[i] * r[i];
        if (sqrt(rs_new) < solver_config.residual_limit)
        {
            printf("* CG solver converged in %d iterations.\n", k + 1);
            break;
        }
        // Update Gnuplot with the current iteration and residual
        update_plot(gnuplot, k, rs_new);
        double beta = rs_new / rs_old;
        for (int i = 0; i < n; i++)
            p[i] = r[i] + beta * p[i];
        rs_old = rs_new;
    }
    if (k == solver_config.max_iteration)
        printf("! CG solver reach to the max of iteration (%d).\n", solver_config.max_iteration);
    // Finalize plotting
    finalize_plot(gnuplot);
    free(r);
    free(p);
    free(Ap);
}


