#include "mkl_iterative_solver.h"

void iterative_preprocess(struct IterativeSolver *solver, const int n, const int *row_ptr, const int *col_idx)
{
    printf("preprocess\n");
    solver->n = n;
    solver->row_ptr = row_ptr;
    solver->col_idx = col_idx;
}

void iterative_analyze(struct IterativeSolver *solver, const int n, const double *val)
{
    printf("analyse\n");
    solver->val = val;
}

void dcg_solve(struct IterativeSolver *solver, const int n, double *x, const double *b)
{
    MKL_INT rci_request, itercount, ipar[128];
    double dpar[128];
    double *tmp = (double *)malloc(sizeof(double) * 4 * n);
    double *temp = (double *)malloc(sizeof(double) * n);
    char matdes[3] = {'d', 'l', 'n'};
    struct matrix_descr descrA;
    sparse_matrix_t csrA;
    sparse_operation_t transA = SPARSE_OPERATION_NON_TRANSPOSE;
    // create sparse matrix
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, n, n, const_cast<int *>(solver->row_ptr), const_cast<int *>(solver->row_ptr) + 1, const_cast<int *>(solver->col_idx), const_cast<double *>(solver->val));

    // init rci
    dcg_init(&n, x, b, &rci_request, ipar, dpar, tmp);

    // set parameters
    ipar[4] = 1000;
    ipar[7] = 0;
    ipar[8] = 1;
    ipar[9] = 1;
    ipar[10] = 0;
    dpar[0] = 1e-10;

    dcg_check(&n, x, b, &rci_request, ipar, dpar, tmp);
    if (rci_request != 0 && rci_request != -1001)
    {
        printf("RCI error: %d\n", rci_request);
    }
    double euclidean_norm;
    // iterate
    while (1)
    {
        dcg(&n, x, b, &rci_request, ipar, dpar, tmp);

        if (rci_request == 0)
        {
            // succeed
            printf("succeed.\n");
            break;
        }
        else if (rci_request == 1)
        {
            // calculate the result
            mkl_sparse_d_mv(transA, 1.0, csrA, descrA, tmp, 0.0, tmp + n);
        }
        else if (rci_request == 2)
        {
            double eone = -1.E0;
            MKL_INT ione = 1;
            double b_norm = dnrm2(&n, b, &ione);
            double loss = dpar[4] / b_norm;
            printf("norm: %e\n", loss);
            if (loss < 1e-10)
                break;
        }
        else if (rci_request == 3)
        {
        }
        else
        {
            // 处理错误
            printf("RCI error: %d\n", rci_request);
            break;
        }
    }
}

void dfgmres_solve(struct IterativeSolver *solver, const int n, double *x, const double *b)
{
    MKL_INT rci_request, itercount, ipar[128];
    double dpar[128];
    double *tmp = (double *)malloc(sizeof(double) * (n * (2 * 150 + 1) + (150 * (150 + 9)) / 2 + 1));
    double *temp = (double *)malloc(sizeof(double) * n);
    char matdes[3] = {'d', 'l', 'n'};
    double *residual = (double *)malloc(sizeof(double) * n);
    double dvar;
    struct matrix_descr descrA;
    sparse_matrix_t csrA;
    sparse_operation_t transA = SPARSE_OPERATION_NON_TRANSPOSE;
    // create sparse matrix
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    for (int i = 0; i < 128; i++)
        ipar[i] = 0;

    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, n, n, const_cast<int *>(solver->row_ptr), const_cast<int *>(solver->row_ptr) + 1, const_cast<int *>(solver->col_idx), const_cast<double *>(solver->val));
    // init rci
    dfgmres_init(&n, x, b, &rci_request, ipar, dpar, tmp);

    // set parameters
    ipar[8] = 1;
    ipar[9] = 0;
    ipar[11] = 1;
    dpar[0] = 1.0E-5;
    dfgmres_check(&n, x, b, &rci_request, ipar, dpar, tmp);
    if (rci_request != 0 && rci_request != -1001)
    {
        printf("RCI error: %d\n", rci_request);
    }
    // iterate
    while (1)
    {
        dfgmres(&n, x, const_cast<double *>(b), &rci_request, ipar, dpar, tmp);

        if (rci_request == 0)
        {
            // succeed
            printf("succeed.\n");
            break;
        }
        else if (rci_request == 1)
        {
            // calculate the result
            mkl_sparse_d_mv(transA, 1.0, csrA, descrA, &tmp[ipar[21] - 1], 0.0, &tmp[ipar[22] - 1]);
            // } else if (rci_request == 2) {
            //     dfgmres_get(&n, x, temp, &rci_request, ipar, dpar, tmp, &itercount);
            //     mkl_sparse_d_mv(transA, 1.0, csrA, descrA, temp, 0.0, residual);
            //     dvar  = -1.0E0;
            //     int i = 1;
            //     daxpy(&n, &dvar, b, &i, residual, &i);
            //     dvar = dnrm2(&n, residual, &i);
            //     printf("norm: %e\n", dvar);
            //     if (dvar < 1e-3) break;
            // } else if (rci_request == 3) {
            // } else if (rci_request == 4) {
            //     if (dpar[6] < 1.0E-12) break;
        }
        else
        {
            // 处理错误
            printf("RCI error: %d\n", rci_request);
            break;
        }
    }
    dfgmres_get(&n, x, const_cast<double *>(b), &rci_request, ipar, dpar, tmp, &itercount);
}

void iterative_solve(struct IterativeSolver *solver, const int n, const double *x, const double *b) { dcg_solve(solver, n, x, b); }
