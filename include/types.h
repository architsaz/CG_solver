#ifndef TYPE_H
#define TYPE_H
#include <stdbool.h>
 
typedef struct {
	int max_iteration;
	double residual_limit;
    bool showplot;
}SolverConfig;

// Sparse Matrix in CSR Format
typedef struct
{
    int *row_ptr;
    int *col_index;
    double *values;
    int n;   // Number of rows (or columns for square matrix)
    int nnz; // Number of non-zero values
} CSRMatrix;


#endif
