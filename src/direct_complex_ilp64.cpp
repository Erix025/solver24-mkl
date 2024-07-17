#include "mkl_direct_solver.h"
#include <stdio.h>
//! complex system
void direct_preprocess_complex_ilp64(struct DirectComplexSolver_ILP64 *solver, const long long int n, const long long int *row_ptr, const long long int *col_idx)
{
    solver->n = n;
    solver->row_ptr = row_ptr;
    solver->col_idx = col_idx;
    solver->mtype = 13; // matrix type
    long long int maxfct, mnum, phase, error;

    for (int i = 0; i < 64; i++)
    {
        solver->pt[i] = 0;
        solver->iparm[i] = 0;
    }
    solver->iparm[0] = 1;
    solver->iparm[1] = 2;
    solver->iparm[3] = 91;
    solver->iparm[7] = 10;
    solver->iparm[9] = 20;
    solver->iparm[20] = 0;
    solver->iparm[23] = 0;
    solver->iparm[34] = 1;
    solver->iparm[59] = 1;

    maxfct = 1;
    mnum = 1;
    solver->phase = 11; // analysis phase
    error = 0;
    long long int nrhs = 1;
    long long int msglvl = 1;

    printf("Preprocess start.\n");
    // initial
    pardiso_64(solver->pt, &maxfct, &mnum, &solver->mtype, &solver->phase, &n, NULL, row_ptr, col_idx, NULL, &nrhs, solver->iparm, &msglvl, NULL, NULL,
               &error);
    if (error != 0)
    {
        printf("PARDISO error: %d\n", error);
    }
    printf("Preprocess End.\n");
}

void direct_analyze_complex_ilp64(struct DirectComplexSolver_ILP64 *solver, const long long int n, const double *val, const double *val_im)
{
    long long int maxfct, mnum, phase, error;

    maxfct = 1;
    mnum = 1;
    error = 0;
    long long int nrhs = 1;
    long long int msglvl = 1;

    // construct complex array
    size_t nnz = solver->row_ptr[solver->n];
    solver->val = (MKL_Complex16 *)malloc(sizeof(MKL_Complex16) * nnz);
    for (size_t i = 0; i < nnz; i++)
    {
        solver->val[i].real = val[i];
        solver->val[i].imag = val_im[i];
    }

    solver->phase = 22; // analysis phase
    printf("Analysis process start.\n");
    pardiso_64(solver->pt, &maxfct, &mnum, &solver->mtype, &solver->phase, &n, (void *)solver->val, solver->row_ptr, solver->col_idx, NULL, &nrhs,
               solver->iparm, &msglvl, NULL, NULL, &error);
    printf("Analysis process end.\n");

    if (error != 0)
    {
        printf("PARDISO error: %d\n", error);
    }
}

void direct_solve_complex_ilp64(struct DirectComplexSolver_ILP64 *solver, const long long int n, double *x, double *x_im, const double *b, const double *b_im)
{
    long long int maxfct, mnum, phase, error;

    maxfct = 1;
    mnum = 1;
    error = 0;
    long long int msglvl = 1;
    long long int nrhs = 1;

    // construct x and b
    solver->x = (MKL_Complex16 *)malloc(sizeof(MKL_Complex16) * n);
    solver->b = (MKL_Complex16 *)malloc(sizeof(MKL_Complex16) * n);
    for (size_t i = 0; i < n; i++)
    {
        solver->x[i].real = x[i];
        solver->x[i].imag = x_im[i];
        solver->b[i].real = b[i];
        solver->b[i].imag = b_im[i];
    }

    solver->phase = 33; // analysis phase
    printf("Solve process start.\n");

    pardiso_64(solver->pt, &maxfct, &mnum, &solver->mtype, &solver->phase, &n, solver->val, solver->row_ptr, solver->col_idx, NULL, &nrhs, solver->iparm,
               &msglvl, solver->b, (double *)solver->x, &error);
    printf("Solve process End.\n");

    if (error != 0)
    {
        printf("PARDISO error: %d\n", error);
    }

    // unpack mkl complex

    for (size_t i = 0; i < n; i++)
    {
        x[i] = solver->x[i].real;
        x_im[i] = solver->x[i].imag;
    }
}
void direct_release_complex_ilp64(struct DirectComplexSolver_ILP64 *solver)
{
    printf("Release.\n");
    long long int maxfct, mnum, error;
    maxfct = 1;
    mnum = 1;
    error = 0;
    long long int msglvl = 1;
    long long int nrhs = 1;
    // release resources
    solver->phase = -1;
    pardiso_64(solver->pt, &maxfct, &mnum, &solver->mtype, &solver->phase, &solver->n, solver->val, solver->row_ptr, solver->col_idx, NULL, &nrhs, solver->iparm,
               &msglvl, NULL, NULL, &error);
    if (error != 0)
    {
        printf("PARDISO error: %d\n", error);
    }
}