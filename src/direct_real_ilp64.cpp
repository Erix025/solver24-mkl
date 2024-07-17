#include "mkl_direct_solver.h"
#include <stdio.h>
//! real system
void direct_preprocess_ilp64(struct DirectSolver_ILP64 *solver, const long long int n, const long long int *row_ptr, const long long int *col_idx)
{
    solver->n = n;
    solver->row_ptr = row_ptr;
    solver->col_idx = col_idx;
    solver->mtype = 11; // matrix type
    long long int maxfct, mnum, phase, error;
    for (int i = 0; i < 64; i++)
    {
        solver->pt[i] = 0;
        solver->iparm[i] = 0;
    }
    solver->iparm[0] = 1;
    solver->iparm[1] = 3;
//    solver->iparm[3] = 91;
//    solver->iparm[7] = 10;
//    solver->iparm[9] = 20;
    solver->iparm[20] = 0;
    solver->iparm[23] = 0;
    solver->iparm[34] = 1;
    solver->iparm[59] = 1;
    if (solver->enable_ordering)
        solver->iparm[4] = 1;
    else
        solver->iparm[4] = 2;

    maxfct = 1;
    mnum = 1;
    solver->phase = 11; // analysis phase
    error = 0;
    long long int nrhs = 1;
    long long int msglvl = 1;

    printf("Preprocess start.\n");
    // initial
    pardiso_64(solver->pt, &maxfct, &mnum, &solver->mtype, &solver->phase, &solver->n, NULL, row_ptr, col_idx, solver->parm, &nrhs, solver->iparm, &msglvl, NULL, NULL,
               &error);
    if (error != 0)
    {
        printf("PARDISO error: %d\n", error);
    }
    printf("Preprocess End.\n");
}

void direct_analyze_ilp64(struct DirectSolver_ILP64 *solver, const long long int n, const double *val)
{
    solver->val = val;
    long long int maxfct, mnum, phase, error;

    maxfct = 1;
    mnum = 1;
    error = 0;
    long long int nrhs = 1;
    long long int msglvl = 1;

    solver->phase = 22; // analysis phase
    printf("Analysis process start.\n");
    pardiso_64(solver->pt, &maxfct, &mnum, &solver->mtype, &solver->phase, &solver->n, (void *)solver->val, solver->row_ptr, solver->col_idx, NULL, &nrhs,
               solver->iparm, &msglvl, NULL, NULL, &error);
    printf("Analysis process end.\n");

    if (error != 0)
    {
        printf("PARDISO error: %d\n", error);
    }
}

void direct_solve_ilp64(struct DirectSolver_ILP64 *solver, const long long int n, const double *x, const double *b)
{
    long long int maxfct, mnum, error;

    maxfct = 1;
    mnum = 1;
    error = 0;
    long long int msglvl = 1;
    long long int nrhs = solver->num_b;

    solver->phase = 33; // analysis phase
    printf("Solve process start.\n");
    pardiso_64(solver->pt, &maxfct, &mnum, &solver->mtype, &solver->phase, &solver->n, solver->val, solver->row_ptr, solver->col_idx, NULL, &nrhs, solver->iparm, &msglvl, (double *)b, (double *)x, &error);
    printf("Solve process End.\n");

    if (error != 0)
    {
        printf("PARDISO error: %d\n", error);
    }
}
void direct_release_ilp64(struct DirectSolver_ILP64 *solver)
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
    pardiso_64(solver->pt, &maxfct, &mnum, &solver->mtype, &solver->phase, &solver->n, solver->val, solver->row_ptr, solver->col_idx, NULL, &nrhs, solver->iparm, &msglvl, NULL, NULL, &error);
    if (error != 0)
    {
        printf("PARDISO error: %d\n", error);
    }
}
