#include <stdio.h>
#include <stdlib.h>
#include "CGSolver.h"

void generate_esure(int grid_width, int grid_height, int *esure) {
    int num_elements = grid_width * grid_height;
    for (int i = 0; i < num_elements; i++) {
        int row = i / grid_width;
        int col = i % grid_width;

        // Determine neighbors
        int bottom = (row > 0) ? (i - grid_width) : -1;  // Bottom boundary (-1)
        int right = (col < grid_width - 1) ? (i + 1) : -2;  // Right boundary (-2)
        int top = (row < grid_height - 1) ? (i + grid_width) : -3;  // Top boundary (-3)
        int left = (col > 0) ? (i - 1) : -4;  // Left boundary (-4)

        // Store neighbors in the esure array
        esure[i * 4 + 0] = bottom;
        esure[i * 4 + 1] = right;
        esure[i * 4 + 2] = top;
        esure[i * 4 + 3] = left;
    }
}

int main()
{
    // mesh
    int grid_width = 1000;
    int grid_height = 1000;
    int nelem = grid_height * grid_width;
    //int npoin = 25;
    int nrpoin = 4;
    int *esure = (int *)malloc((size_t)nelem * 4 * sizeof(int));
    if (esure == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }
    generate_esure(grid_width, grid_height, esure);
    #ifdef DEBUG
        // Print the esure array
        printf("int esure[] = {\n");
        for (int i = 0; i < nelem; i++) {
            printf("    %d, %d, %d, %d", esure[i * 4 + 0], esure[i * 4 + 1], esure[i * 4 + 2], esure[i * 4 + 3]);
            if (i < nelem - 1) {
                printf(",\n");
            }
        }
        printf("\n};\n");
    #endif
    // define sparse matrix of coefficient 
    int *row_ptr = (int *)calloc((size_t)(nelem + 1), sizeof(int));
    int max_nnz = nelem * 6; // At most 4 non-zero entries per row
    double *val = (double *)malloc((size_t)max_nnz * sizeof(double));
    int *col_ind = (int *)malloc((size_t)max_nnz * sizeof(int));
    int nnz = 0;
//    double *coeff = calloc((size_t)nelem * (size_t)nelem, sizeof(double));
    double *RHS = calloc((size_t)nelem, sizeof(double));

    for (int ele = 0; ele < nelem; ele++) {
        nnz++;
        int IDele = nnz-1;
        for (int nei = 0; nei < nrpoin; nei++) {
            int neighbor = esure[nrpoin * ele + nei];
            if (neighbor >= 0) { // Internal connection
                //coeff[nelem * ele + ele] += 1;               // Subtract full contribution for diagonal
                val[IDele] += 1;
                col_ind[IDele] = ele; 
                //coeff[nelem * ele + neighbor] -= 1;          // Add contribution to the neighbor
                val[nnz] = -1;
                col_ind[nnz]= neighbor;
                nnz++;
            } else if (neighbor == -1)
            { // Boundary condition
               // coeff[nelem * ele + ele] += 2;               // Subtract halved contribution (distance halved)
                val[IDele] += 2;
                RHS[ele] += 5 * 2;                          // Add scaled boundary condition to RHS
            } else{
            { // Boundary condition
               // coeff[nelem * ele + ele] += 2;               // Subtract halved contribution (distance halved)
                val[IDele] += 2;
                RHS[ele] += 1 * 2;                          // Add scaled boundary condition to RHS
            }
            }
        }
        row_ptr[ele+1]=nnz;
    }
    row_ptr[nelem]=nnz;
    
    // Conver to CSR format
    #ifdef DEBUG
        //Print the CRS representation
        printf("Values (val): ");
        for (int i = 0; i < nnz; i++) {
            printf("%.2lf ", val[i]);
        }
        printf("\n");
        printf("Column Indices (col_ind): ");
        for (int i = 0; i < nnz; i++) {
            printf("%d ", col_ind[i]);
        }
        printf("\n");
        printf("Row Pointers (row_ptr): ");
        for (int i = 0; i <= nelem; i++) {
            printf("%d ", row_ptr[i]);
        }
        printf("\n");
    #endif
    // Structure matrix A (Coefficient matrix) in CSR format 
    CSRMatrix A;
    A.n = nelem;
    A.nnz = nnz;
    A.row_ptr = row_ptr;
    A.col_index = col_ind;
    A.values = val;

    // unknown vector 
    double *u = (double *)calloc((size_t)nelem,sizeof(double));

    // Solve using CG
    SolverConfig config = {100000,1e-8,true};
    solver_set_config(config); 
    //precond_conjugate_gradient(&A, RHS, u);
    conjugate_gradient(&A, RHS, u);
    //   Output solution
    // printf("Solution u:\n");
    // for (int i = 0; i < grid_height; i++){
    //     for (int j = 0; j < grid_width; j++)
    //     printf("%f ", u[grid_height*i+j]);
    //     printf("\n");
    // }

    free(val);
    free(col_ind);
    free(row_ptr);
    free (RHS);
    free(u);
    free(esure);
    return 0;
}

