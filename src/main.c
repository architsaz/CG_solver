#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "CGSolver.h"
#include "CRSmatfuncs.h"
#include "CDSolver.h"

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
void generate_cell_stat(int nelem,int *esure,int nrpoin,int *cell_stat) {
    for (int ele=0;ele<nelem;ele++){
        for (int i=0;i<nrpoin;i++){
            int neighbor = esure [nrpoin*ele+i];
            if (neighbor<0){
                cell_stat [ele] =  neighbor;
                break;
            }
            
        }
    }
}
int main()
{
    // clock 
    clock_t start_time, end_time;
    double cpu_time_used;
    // mesh
    int grid_width = 3;
    int grid_height = 3;
    int nelem = grid_height * grid_width;
    //int npoin = 25;
    int nrpoin = 4;
    int *esure = (int *)malloc((size_t)nelem * 4 * sizeof(int));
    int *cell_stat = calloc((size_t)nelem,sizeof(int));
    if (esure == NULL || cell_stat == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }
    generate_esure(grid_width, grid_height, esure);
    generate_cell_stat(nelem,esure,nrpoin,cell_stat);
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
        printf("int cell_stat:\n");
        for (int i = 0; i < grid_height; i++){
            for (int j = 0; j < grid_width; j++)
            printf("%d ", cell_stat[grid_height*i+j]);
            printf("\n");
        }
    #endif
    // define sparse matrix of coefficient 
    int *row_ptr = (int *)calloc((size_t)(nelem + 1), sizeof(int));
    int max_nnz = nelem * 6; // At most 4 non-zero entries per row
    double *val = (double *)malloc((size_t)max_nnz * sizeof(double));
    int *col_ind = (int *)malloc((size_t)max_nnz * sizeof(int));
    int nnz = 0;
    double *RHS = (double *)calloc((size_t)nelem, sizeof(double));
    if (row_ptr == NULL || val == NULL || col_ind == NULL || RHS == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }
    // define coefficient matrix in CRS format and right-hand-side vecotr for the system of equation (Au=RHS) 
    
    // for (int ele = 0; ele < nelem; ele++) {
    // nnz++;
    // int IDele = nnz-1;
    //     if (cell_stat[ele] == 0) { // Internal connection
    //         //coeff[nelem * ele + ele] += 1;               // Subtract full contribution for diagonal
    //         val[IDele] = nrpoin;
    //         col_ind[IDele] = ele; 
    //         //coeff[nelem * ele + neighbor] -= 1;          // Add contribution to the neighbor
    //         for (int nei = 0; nei < nrpoin; nei++) {
    //             val[nnz] = -1;
    //             col_ind[nnz]= esure[nrpoin * ele + nei];
    //             nnz++;
    //         }
    //     } else if (cell_stat[ele] == -1)
    //     { // Boundary condition
    //         // coeff[nelem * ele + ele] += 2;               // Subtract halved contribution (distance halved)
    //         val[IDele] = 1;
    //         col_ind[IDele] = ele;
    //         RHS[ele] = 5 ;                          // Add scaled boundary condition to RHS
    //     } else{
    //         // Boundary condition
    //         // coeff[nelem * ele + ele] += 2;               // Subtract halved contribution (distance halved)
    //         val[IDele] = 1;
    //         col_ind[IDele] = ele;
    //         RHS[ele] = 1 ;                          // Add scaled boundary condition to RHS
    //     }
    // row_ptr[ele+1]=nnz;
    // }
    // row_ptr[nelem]=nnz;
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
        printf("RHS: ");
        for (int i = 0; i < nelem; i++) {
            printf("%lf ", RHS[i]);
        }
        printf("\n");
    #endif
    //Check the Positive Difinte conditions for a symetric off-diagonal sparce matrix  
    if (isPositiveDefinite(nelem, val, col_ind, row_ptr)) {
        printf("The matrix is positive definite.\n");
    } else {
        fprintf(stderr,"The Coefficient matrix is not positive definite.\n");
        //return 1;
    }
    // Structure matrix A (Coefficient matrix) in CSR format 
    // int Tn = 3;
    // double Tvalues[] = {4, 1, 1, 1, 3, 1, 2};
    // int Tcolumns[] = {0, 1, 2, 0, 1, 0, 2};
    // int Trow_ptr[] = {0, 3, 5, 7};
    // double Tb[] = {7, 4, 5};
    // int Tnnz = 7;

    // int Tn = 4;
    // double Tvalues[] = {4, 1, 1, 3, 1, 1, 2, 1, 1, 2};
    // int Tcolumns[] =   {0, 1, 0, 1, 2, 1, 2, 3,  2, 3};
    // int Trow_ptr[] = {0, 2, 5, 8, 10};
    // double Tb[] = {6,10,12,11};
    // int Tnnz = 10;

    // int Tn = 5;
    // double Tvalues[] = {4, 1, 1, 3, 1, 1, 2, 1, 1, 3, 1, 1, 2};
    // int Tcolumns[] =   {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4};
    // int Trow_ptr[] = {0, 2, 5, 8, 11, 13};
    // double Tb[] = {6,10,12,20,14};
    // int Tnnz = 13;

    int Tn = 6;
    double Tvalues[] = {6, 2, 2, 5, 2, 2, 4, 1, 1, 4, 2, 2, 5, 3, 3, 6};
    int Tcolumns[] =   {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5};
    int Trow_ptr[] = {0, 2, 5, 8, 11, 14, 16};
    double Tb[] = {8, 9, 7, 7, 10, 9};
    int Tnnz = 16;


    CRSMatrix A;
    A.n = Tn;
    A.nnz = Tnnz;
    A.row_ptr = Trow_ptr;
    A.col_index = Tcolumns;
    A.values = Tvalues;

    // CRSMatrix A;
    // A.n = nelem;
    // A.nnz = nnz;
    // A.row_ptr = row_ptr;
    // A.col_index = col_ind;
    // A.values = val;

    // unknown vector 
    double *u = (double *)calloc((size_t)A.n,sizeof(double));

    // // Solve using CG
    // SolverConfig config = {100000,1e-8,true};
    // solver_set_config(config); 
    // //precond_conjugate_gradient(&A, RHS, u);
    // start_time = clock();
    // conjugate_gradient(&A, Tb, u);
    // end_time = clock();
    // cpu_time_used = (double)(end_time-start_time)/CLOCKS_PER_SEC;
    // printf("! CG Solver execution time : %.2f seconds\n",cpu_time_used);

    // Solve using Cholesky Decomposition 
    start_time = clock();
    solve_cholesky(&A, Tb, u);
    end_time = clock();
    cpu_time_used = (double)(end_time-start_time)/CLOCKS_PER_SEC;
    printf("! Cholesky Solver execution time : %.2f seconds\n",cpu_time_used);

    //Output solution
    printf("Solution u:\n");
    for (int i = 0; i < A.n; i++)
        printf("%f \n", u[i]);
    

    //Output solution
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
    free(cell_stat);
    return 0;
}

