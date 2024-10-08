/*
 * This is the reference iterative method solver code of SolverChallenge24
 * author: Li Zhao and Haibing Sun
 */

// system file
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>

// read matrix file
#include "mmio_highlevel.h"

// utlise function file
#include "utlise.h"
#include "utlise_long.h"

// iterative solver code
#include "mkl_direct_solver.h"
#include "mkl_iterative_solver.h"
int main(int argc, char** argv)
{
    int     m, n, nnzA, isSymmetricA;
    int*    row_ptr = NULL;         // the csr row pointer array of matrix A
    int*    col_idx = NULL;         // the csr column index array of matrix A
    double* val     = NULL;         // the csr value array of matrix A (real number)
    double* val_im  = NULL;         // the csr value array of matrix A (imaginary number)
    double *x = NULL, *x_im = NULL; // solution vector x, (x: real number, x_im: imaginary number)
    double *b = NULL, *b_im = NULL; // right-hand side vector b, (b: real number, b_im: imaginary number)
    double  tt, time;

    char* filename_matrix;      // the filename of matrix A
    char* filename_b;           // the filename of right-hand side vector b
    int   read_matrix_base = 1; // 0-base or 1-base, default 1-base
    int   type             = 0; // type to output time, 0: end to end time; 1:solver time + solve time; 2:solve time; default 0
    int   test_frequency   = 1; // run code frequency
    int   sys_type         = 0; // type of algebraic systems, 0: real, 1: complex; default 0

    /* ========================================== */
    // Step 0: Read command line argument
    /* ========================================== */
    if (argc == 3) {
        filename_matrix = argv[1];
        filename_b      = argv[2];
    } else if (argc == 4) {
        filename_matrix  = argv[1];
        filename_b       = argv[2];
        read_matrix_base = atoi(argv[3]);
        if (!((read_matrix_base == 0) || (read_matrix_base == 1))) {
            fprintf(stdout, "read_base error input %d\n", read_matrix_base);
            print_help();
            return -1;
        }
    } else if (argc == 5) {
        filename_matrix  = argv[1];
        filename_b       = argv[2];
        read_matrix_base = atoi(argv[3]);
        if (!((read_matrix_base == 0) || (read_matrix_base == 1))) {
            fprintf(stdout, "read_base error input %d\n", read_matrix_base);
            print_help();
            return -1;
        }
        type = atoi(argv[4]);
        if (type > 2 || type < 0) {
            fprintf(stdout, "type error input %d\n", type);
            print_help();
            return -1;
        }
    } else if (argc == 6) {
        filename_matrix  = argv[1];
        filename_b       = argv[2];
        read_matrix_base = atoi(argv[3]);
        if (!((read_matrix_base == 0) || (read_matrix_base == 1))) {
            fprintf(stdout, "read_base error input %d\n", read_matrix_base);
            print_help();
            return -1;
        }
        type = atoi(argv[4]);
        if (type > 2 || type < 0) {
            fprintf(stdout, "type error input %d\n", type);
            print_help();
            return -1;
        }
        sys_type = atoi(argv[5]);
        if (sys_type > 1 || sys_type < 0) {
            fprintf(stdout, "sys_type error input %d\n", sys_type);
            print_help();
            return -1;
        }
    } else {
        print_help();
        return -1;
    }

    //! If the user does not specify parameter sys_type,
    //! it will automatically determine whether it is a real system or a complex system according to the mtx file
    if (argc < 6) read_mtx_header(filename_matrix, &sys_type);

    fprintf(stdout, "matrix name :      %s\nvectorb name :     %s\nread function :    base-%d\ntype :             %d\nsys_type :         %d\n",
            filename_matrix, filename_b, read_matrix_base, type, sys_type);

    /* ========================================== */
    // Step 1: Load matrix and rhs from mtx files
    /* ========================================== */
    if (sys_type == 0) // real system
    {
        // load matrix
        mmio_allinone(&m, &n, &nnzA, &isSymmetricA, &read_matrix_base, &row_ptr, &col_idx, &val, filename_matrix);
        if (m != n) {
            fprintf(stdout, "Invalid matrix size.\n");
            return 0;
        }

        x = (double*)malloc(sizeof(double) * n);
        b = (double*)malloc(sizeof(double) * n);

        // load right-hand side vector b
        load_vector(n, b, filename_b);

    } else { // complex system
        mmio_allinone_complex(&m, &n, &nnzA, &isSymmetricA, &read_matrix_base, &row_ptr, &col_idx, &val, &val_im, filename_matrix);
        if (m != n) {
            fprintf(stdout, "Invalid matrix size.\n");
            return 0;
        }

        x    = (double*)malloc(sizeof(double) * n);
        x_im = (double*)malloc(sizeof(double) * n);
        b    = (double*)malloc(sizeof(double) * n);
        b_im = (double*)malloc(sizeof(double) * n);

        // load right-hand side vector b
        load_b_complex(n, b, b_im, filename_b);
    }

    /* ========================================== */
    // Step 2: Solve the linear system
    /* ========================================== */
    if (sys_type == 0) // real system
    {
        // iterative solver sample
        struct IterativeSolver iterative_solver;
        struct DirectSolver    direct_solver;
        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++) {
                // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                memset(x, 0.0, sizeof(double) * n); // initial vector x
                direct_preprocess(&direct_solver, n, row_ptr, col_idx);
                direct_analyze(&direct_solver, n, val);
                direct_solve(&direct_solver, n, x, b);
                iterative_preprocess(&iterative_solver, n, row_ptr, col_idx);
                iterative_analyze(&iterative_solver, n, val);
                iterative_solve(&iterative_solver, n, x, b);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        } else if (type == 1) // check preprocess + iterative_solver time
        {
            direct_preprocess(&direct_solver, n, row_ptr, col_idx);
            iterative_preprocess(&iterative_solver, n, row_ptr, col_idx);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++) {
                // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                memset(x, 0.0, sizeof(double) * n); // initial vector x

                direct_analyze(&direct_solver, n, val);
                direct_solve(&direct_solver, n, x, b);
                iterative_analyze(&iterative_solver, n, val);
                iterative_solve(&iterative_solver, n, x, b);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        } else { // check iterative_solver time

            direct_preprocess(&direct_solver, n, row_ptr, col_idx);
            direct_analyze(&direct_solver, n, val);
            iterative_preprocess(&iterative_solver, n, row_ptr, col_idx);
            iterative_analyze(&iterative_solver, n, val);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++) {
                // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                memset(x, 0.0, sizeof(double) * n); // initial vector x
                direct_solve(&direct_solver, n, x, b);
                iterative_solve(&iterative_solver, n, x, b);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }

    } else { // complex system
        // iterative solver sample
        struct IterativeComplexSolver mysolver;
        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++) {
                // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                memset(x, 0.0, sizeof(double) * n); // initial vector x
                memset(x_im, 0.0, sizeof(double) * n);

                iterative_preprocess_complex(&mysolver, n, row_ptr, col_idx);
                iterative_analyze_complex(&mysolver, n, val, val_im);
                iterative_solve_complex(&mysolver, n, x, x_im, b, b_im);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        } else if (type == 1) // check preprocess + iterative_solver time
        {
            iterative_preprocess_complex(&mysolver, n, row_ptr, col_idx);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++) {
                // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                memset(x, 0.0, sizeof(double) * n); // initial vector x
                memset(x_im, 0.0, sizeof(double) * n);

                iterative_analyze_complex(&mysolver, n, val, val_im);
                iterative_solve_complex(&mysolver, n, x, x_im, b, b_im);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        } else { // check iterative_solver time

            iterative_preprocess_complex(&mysolver, n, row_ptr, col_idx);
            iterative_analyze_complex(&mysolver, n, val, val_im);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++) {
                // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                memset(x, 0.0, sizeof(double) * n); // initial vector x
                memset(x_im, 0.0, sizeof(double) * n);

                iterative_solve_complex(&mysolver, n, x, x_im, b, b_im);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
    }

    fprintf(stdout, "------------------------------------------\n");

    /* ========================================== */
    // Step 3: Check time, memory and correctness
    /* ========================================== */
    // check time
    if (type == 0) {
        fprintf(stdout, "CHECK end to end time :         %12.6lf ms\n", time);
    } else if (type == 1) {
        fprintf(stdout, "CHECK solver + solve time :     %12.6lf ms\n", time);
    } else if (type == 2) {
        fprintf(stdout, "CHECK solve time :              %12.6lf ms\n", time);
    }

    // check the memory
    mem_usage();

    // store x to a file
    char* answer_x = "answer_x.rhs";
    if (sys_type == 0) // real rhs
        store_x(n, x, answer_x);
    else // complex rhs
        store_x_complex(n, x, x_im, answer_x);

    // check the correctness
    if (sys_type == 0) // real system
    {
        // 1) using double precision
        check_correctness(n, row_ptr, col_idx, val, x, b);
        // 2) using long double precision
        check_correctness_ld_d2ld(n, row_ptr, col_idx, val, x, b);

    } else { // complex system
        // 1) using double precision
        check_correctness_complex(n, row_ptr, col_idx, val, val_im, x, x_im, b, b_im);
        // 2) using long double precision
        check_correctness_complex_ld_d2ld(n, row_ptr, col_idx, val, val_im, x, x_im, b, b_im);

        free(val_im);
        free(x_im);
        free(b_im);
    }

    free(row_ptr);
    free(col_idx);
    free(val);
    free(x);
    free(b);

    return 0;
}
