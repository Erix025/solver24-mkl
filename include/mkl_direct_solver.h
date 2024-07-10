#include "mkl.h"
#undef MKL_ILP64
// please add your code in this file
struct DirectSolver
{
    MKL_INT n; // size of matrix
    // CSR format
    const int *row_ptr;
    const int *col_idx;
    int num_b;
    const double *val;
    void *pt[64];
    MKL_INT iparm[64];
    MKL_INT mtype;
    MKL_INT phase;
};

struct DirectSolver_ILP64
{
    long long int n; // size of matrix
    // CSR format
    const long long int *row_ptr;
    const long long int *col_idx;
    int num_b;
    const double *val;
    void *pt[64];
    long long int iparm[64];
    long long int mtype;
    long long int phase;
};

struct DirectComplexSolver
{
    // add struct array
    int n;
    const int *row_ptr;
    const int *col_idx;
    MKL_Complex16 *val;
    void *pt[64];
    MKL_INT iparm[64];
    MKL_INT mtype;
    MKL_INT phase;
    MKL_Complex16 *x;
    MKL_Complex16 *b;
};

//! real system
void direct_preprocess(struct DirectSolver *solver, const int n, const int *row_ptr, const int *col_idx);

void direct_analyze(struct DirectSolver *solver, const int n, const double *val);

void direct_solve(struct DirectSolver *solver, const int n, const double *x, const double *b);

void direct_release(struct DirectSolver *solver);

void direct_preprocess_ilp64(struct DirectSolver_ILP64 *solver, const long long int n, const long long int *row_ptr, const long long int *col_idx);

void direct_analyze_ilp64(struct DirectSolver_ILP64 *solver, const long long int n, const double *val);

void direct_solve_ilp64(struct DirectSolver_ILP64 *solver, const long long int n, const double *x, const double *b);

void direct_release_ilp64(struct DirectSolver_ILP64 *solver);

//! complex system
void direct_preprocess_complex(struct DirectComplexSolver *solver, const int n, const int *row_ptr, const int *col_idx);

void direct_analyze_complex(struct DirectComplexSolver *solver, const int n, const double *val, const double *val_im);

void direct_solve_complex(struct DirectComplexSolver *solver, const int n, double *x, double *x_im, const double *b, const double *b_im);

void direct_release_complex(struct DirectComplexSolver *solver);