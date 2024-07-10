#include "mkl_iterative_solver.h"
//! complex system
void iterative_preprocess_complex(struct IterativeComplexSolver *solver, const int n, const int *row_ptr, const int *col_idx)
{
    // add the code
}

void iterative_analyze_complex(struct IterativeComplexSolver *solver, const int n, const double *val, const double *val_im)
{
    // add the code
}

void iterative_solve_complex(struct IterativeComplexSolver *solver, const int n, const double *x, const double *x_im, const double *b,
                             const double *b_im)
{
}
