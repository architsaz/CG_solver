#ifndef CRSMAT_TYPE_H
#define CRSMAT_TYPE_H

// Sparse Matrix in CSR Format
typedef struct
{
    int *row_ptr;
    int *col_index;
    double *values;
    int n;   // Number of rows (or columns for square matrix)
    int nnz; // Number of non-zero values
} CRSMatrix;


#endif
