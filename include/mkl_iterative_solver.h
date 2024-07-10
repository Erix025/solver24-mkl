#include "mkl.h"
#include <stdio.h>
struct IterativeSolver
{
    int n; // size of matrix
    // CSR format
    const int *row_ptr;
    const int *col_idx;
    const double *val;
};

struct IterativeComplexSolver
{
    // add struct array
};

void iterative_preprocess(struct IterativeSolver *solver, const int n, const int *row_ptr, const int *col_idx);

void iterative_analyze(struct IterativeSolver *solver, const int n, const double *val);

void dcg_solve(struct IterativeSolver *solver, const int n, const double *x, const double *b);

void dfgmres_solve(struct IterativeSolver *solver, const int n, const double *x, const double *b);

void iterative_solve(struct IterativeSolver *solver, const int n, const double *x, const double *b);

//! complex system
void iterative_preprocess_complex(struct IterativeComplexSolver *solver, const int n, const int *row_ptr, const int *col_idx);

void iterative_analyze_complex(struct IterativeComplexSolver *solver, const int n, const double *val, const double *val_im);

void iterative_solve_complex(struct IterativeComplexSolver *solver, const int n, const double *x, const double *x_im, const double *b, const double *b_im);