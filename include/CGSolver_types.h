#ifndef CGSOLVER_TYPES_H
#define CGSOLVER_TYPES_H
#include <stdbool.h>

// all configuration parameters for CG Solver 
typedef struct {
	int max_iteration;
	double residual_limit;
    bool showplot;
}SolverConfig;

#endif