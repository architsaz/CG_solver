#include "CRSMat_types.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// print the RRS matrix in terminal
void print_CRSMatrix(CRSMatrix *mat) {
    printf("CRS Matrix (n = %d):\n", mat->n);
    printf("row_ptr: [");
    for (int i = 0; i <= mat->n; i++) {
        printf("%d", mat->row_ptr[i]);
        if (i < mat->n) printf(", ");
    }
    printf("]\n");

    printf("col_index: [");
    for (int i = 0; i < mat->row_ptr[mat->n]; i++) {
        printf("%d", mat->col_index[i]);
        if (i < mat->row_ptr[mat->n] - 1) printf(", ");
    }
    printf("]\n");

    printf("values: [");
    for (int i = 0; i < mat->row_ptr[mat->n]; i++) {
        printf("%.6f", mat->values[i]);
        if (i < mat->row_ptr[mat->n] - 1) printf(", ");
    }
    printf("]\n");
}
// Matrix-Vector Multiplication
void csr_matvec(CRSMatrix *A, double *x, double *y)
{
    for (int i = 0; i < A->n; i++)
    {
        y[i] = 0.0;
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            y[i] += A->values[j] * x[A->col_index[j]];
        }
    }
}
// Function to perform Cholesky decomposition
void choleskyDecomposition(CRSMatrix *A, CRSMatrix *L) {
    // Check for valid input
    if (A == NULL || L == NULL || A->n <= 0) {
        fprintf(stderr, "Invalid input matrices.\n");
        exit(EXIT_FAILURE);
    }
    // Initialize L
    L->n = A->n;
    L->nnz = A->n*(A->n+1)/2;
    L->col_index = (int *)calloc((size_t)L->nnz,sizeof(int));
    L->values = (double *)calloc((size_t)L->nnz,sizeof(double));
    L->row_ptr = (int *)calloc((size_t)(L->n+1),sizeof(int));
    int nnz1=0;
    for (int i = 0;i<L->n;i++) {
        L->row_ptr [i]=nnz1;
        for (int j = 0;j<=i;j++){
            L->col_index [nnz1]=j;
            nnz1++;
        }
    }
    L->row_ptr [L->n]=nnz1;

    for (int j = 0; j < A->n; j++) { // loop over the columns 
        for (int i = j; i < A->n ; i++) { // loop over row
            int n = -1;
            for (int ii = A->row_ptr [i];ii<A->row_ptr [i+1];ii++){
                if (A->col_index[ii]==j) {
                n= ii-A->row_ptr [i];
                break;
                }
            }
            int max_nonzero_in_row = A->row_ptr [i+1] - A->row_ptr [i];
        // for (int n = 0; n < max_nonzero_in_row ; n++) { // loop over columns
        //     int j = A->col_index[A->row_ptr [i]+n]; // column index for each nonzero element
            if (j > i) continue; // make sure picked element in the lower triangle half of Matrix
            double sum = 0.0;

            if (j == i) {  // Diagonal elements
                for (int kk = 0; kk < max_nonzero_in_row; kk++)
                {                   
                    int k = A->col_index[A->row_ptr [i]+kk];
                    if (k >= j) continue;
                    double Lnjk = L->values [L->row_ptr [j]+k];
                    sum += Lnjk * Lnjk;
                }
                double Anjj = (n!=-1) ? A->values [A->row_ptr[j]+n] : 0 ;
                if (Anjj - sum <= 0) {
                    fprintf(stderr,"Matrix is not positive definite (Ajj-sum <= 0).\n");
                    exit(EXIT_FAILURE);
                }
                //double Lnjj = L->values [L->row_ptr[j]+j];
                L->values [L->row_ptr[j]+j] = sqrt(Anjj - sum);

            } else {  // Off-diagonal elements
                for (int kk = 0; kk < max_nonzero_in_row; kk++){
                    int k = A->col_index[A->row_ptr [i]+kk];
                    if (k >= j) continue;    
                    double Lnik = L->values [L->row_ptr [i]+k];
                    double Lnjk = L->values [L->row_ptr [j]+k];            
                    sum += Lnik * Lnjk;
                }
                double Lnjj = L->values [L->row_ptr[j]+j];
                double Anij = (n!=-1) ? A->values [A->row_ptr[i]+n] : 0;
                L->values [L->row_ptr[i]+j] = (Anij - sum) / Lnjj;
            }
        }
    }
}
// Function to perform forward substitution
void forward_substitution(CRSMatrix *L, double *b, double *y) {
    int n = L->n;
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        int diag_col_indx = -1;
        for (int k = L->row_ptr[i]; k < L->row_ptr[i + 1]; k++) {
            int j = L->col_index[k];
            if (j < i) {
                y[i] -= L->values[k] * y[j];
            }
            if (j == i ){
                diag_col_indx = k ;
            }
        }
        if (diag_col_indx == -1){
            fprintf(stderr,"forward substitution can not happen\n* the diagonal element of matrix L at row %d is zero\n",i);
            exit(EXIT_FAILURE);
        }
        //printf("i: %d diag_indx: %d y[i] : %lf lii: %lf\n",i,diag_col_indx,y[i],L->values[diag_col_indx]);
        y[i] /= L->values[diag_col_indx];
    }
}
// Function to perform backward substitution
void backward_substitution(CRSMatrix *L, double *y, double *x) {
    int n = L->n;
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        int diag_col_indx = -1;
        for (int k = L->row_ptr[i]; k < L->row_ptr[i + 1]; k++) {
            int j = L->col_index[k];
            if (j > i) {
                x[i] -= L->values[k] * x[j];
            }
            if (j == i ){
                diag_col_indx = k ;
            }
        }
        if (diag_col_indx == -1){
            fprintf(stderr,"forward substitution can not happen\n* the diagonal element of matrix L at row %d is zero\n",i);
            exit(EXIT_FAILURE);
        }
        x[i] /= L->values[diag_col_indx];
    }
}
// Function for transpose a CRS matrices 
void transCRSmat (CRSMatrix *A, CRSMatrix *AT){
    // Check for valid input
    if (A == NULL || AT == NULL ) {
        fprintf(stderr, "Invalid input matrices.\n");
        exit(EXIT_FAILURE);
    }
    // Initialize the AT CRS Matrix
    AT->nnz = A->nnz;
    AT->n = A->n;
    AT->col_index = (int *)malloc((size_t)AT->nnz*sizeof(int));
    AT->row_ptr = (int *)calloc((size_t)(AT->n+1),sizeof(int));
    AT->values = (double *)malloc((size_t)AT->nnz*sizeof(double));
    
    // Step 1: Count the number of non-zero entries per column in the original matrix
    for (int i = 0; i < A->nnz; i++) {
        AT->row_ptr[A->col_index[i] + 1]++;
    }

    // Step 2: Compute the prefix sum to determine the row_ptr of the transposed matrix
    printf("AT->row_ptr [] : \n%d ",AT->row_ptr[0]);
    for (int i = 1; i <= AT->n; i++) {
        AT->row_ptr[i] += AT->row_ptr[i - 1];
        printf("%d ",AT->row_ptr[i]);
    }
    printf("\n");

    // step 3 : Populate the values and col_indx for the transposed matrix
    int *current_pos = (int *)malloc((size_t)AT->n * sizeof(int));
    for (int i = 0; i < AT->n; i++) {
        current_pos[i] = AT->row_ptr[i];
    }
    printf("k: \n");
    for (int i = 0; i < A->n; i++) {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
            int k = current_pos[A->col_index[j]]; // Get the position for  start of each row in AT
            AT->values[k] = A->values[j];        // Copy the value
            AT->col_index[k] = i;                // The row index of A becomes the column index of AT
            current_pos[A->col_index[j]]++;    // Increment the position for this column
        }
    }
    free(current_pos);

}
// Function to check if the leading k x k submatrix is positive definite using CRS format
int isPositiveDefiniteMinor(int k, double *values, int *columns, int *row_ptr) {
    // Initialize an array to store the diagonal elements of the Cholesky decomposition
    double *L_diag = (double *)malloc((size_t)k * sizeof(double)); // Store only diagonal elements
    if (L_diag == NULL) {
        fprintf(stderr, "Failed to allocate memory for L_diag\n");
        return 0; // Memory allocation failed
    }

    for (int i = 0; i < k; i++) {
        double sum = 0.0;

        // Process the diagonal element
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            int col = columns[j];
            if (col < i) {
                sum -= values[j] * values[j] / L_diag[col]; // Contribution from previous rows
            } else if (col == i) {
                sum += values[j]; // Diagonal element
            }
        }

        // Check if the diagonal element is positive
        if (sum <= 0) {
            free(L_diag);
            fprintf(stderr, "Diagonal element is non-positive at row %d\n", i);
            return 0; // Not positive definite
        }

        L_diag[i] = sqrt(sum); // Update the diagonal element
    }

    free(L_diag);
    return 1; // Positive definite if all diagonal elements are positive
}
// Function to check if a matrix is positive definite based on CRS format
int isPositiveDefinite(int n, double *values, int *columns, int *row_ptr) {
    // Step 1: Check diagonal elements
    for (int i = 0; i < n; i++) {
        int start = row_ptr[i];
        int end = row_ptr[i + 1];
        int found_diagonal = 0;
        for (int j = start; j < end; j++) {
            if (columns[j] == i) {
                if (values[j] <= 0) {
                    fprintf(stderr,"Not positive definite due to diagonal element is %lf <= 0\n",values[j]);
                    return 0; // Not positive definite if diagonal element is <= 0
                }
                found_diagonal = 1;
                break;
            }
        }
        if (!found_diagonal) {
            fprintf(stderr,"Not positive definite due to missing diagonal element at row %d <= 0\n",i);
            return 0; // Not positive definite if diagonal element is missing
        }
    }

    // Step 2: Check leading principal minors using CRS format
    for (int k = 1; k <= n; k++) {
        // Perform Cholesky decomposition directly for leading k x k submatrix
        if (!isPositiveDefiniteMinor(k, values, columns, row_ptr)) {
            fprintf(stderr,"Not positive definite due to leading minor is non-positive\n");
            return 0; // Not positive definite if any leading minor is non-positive
        }
    }

    return 1; // The matrix is positive definite
}
// Function to count non-zero elements in the matrix
int countNonZero(int rows, int cols, double *matrix) {
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (fabs(matrix[cols * i + j]) > 1e-8) {  // Use cols here
                count++;
            }
        }
    }
    return count;
}
// conver sparse Matrix to CRS format 
void convertToCRS(int rows, int cols, double *matrix, double *val, int *col_ind, int *row_ptr) {
    int k = 0;  // Index for `val` and `col_ind`
    row_ptr[0] = 0;  // First row starts at index 0

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (fabs(matrix[cols * i + j]) > 10E-9) {
                val[k] = matrix[cols * i + j];  // Store non-zero value
                col_ind[k] = j;                  // Store column index
                k++;
            }
        }
        row_ptr[i + 1] = k;  // Update row pointer to the total count of non-zero elements so far
    }
}