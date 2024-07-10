/*
 * This is the reference direct method solver code of SolverChallenge24
 * author: Li Zhao and Haibing Sun
 */

// system file
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/resource.h>
#include <sys/time.h>
#include <sstream>
#include <iomanip>
// read matrix file
#include "mmio_highlevel.h"

// utlise function file
#include "utlise.h"
#include "utlise_long.h"

// direct solver code
#include "mkl_direct_solver.h"

struct SolverArgument
{
    std::string filename_matrix;
    std::string filename_b;
    int read_matrix_base = 1; // 0-base or 1-base, default 1-base
    int type = 0;             // type to output time, 0: end to end time; 1:solver time + solve time; 2:solve time; default 0
    int test_frequency = 3;   // run code frequency
    int sys_type = 0;         // type of algebraic systems, 0: real, 1: complex; default 0
    int num_b = 1;            // number of b to solve, default 1
};

struct TimerRec
{
    double preprocess_time = 0;
    double analyze_time = 0;
    double solve_time = 0;
};

using SolverArgument = struct SolverArgument;
using TimerRec = struct TimerRec;

std::string get_index_filename(std::string &prefix, int index)
{
    std::string filename = prefix + "_t" + std::to_string(index) + ".rhs";
    return filename;
}

void parse_args(int argc, char **argv, SolverArgument *args)
{
    if (argc >= 3)
    {
        args->filename_matrix = argv[1];
        args->filename_b = argv[2];
    }
    if (argc >= 4)
    {
        args->read_matrix_base = atoi(argv[3]);
        if (!((args->read_matrix_base == 0) || (args->read_matrix_base == 1)))
        {
            fprintf(stdout, "read_base error input %d\n", args->read_matrix_base);
            print_help();
            exit(-1);
        }
    }
    if (argc >= 5)
    {
        args->type = atoi(argv[4]);
        if (args->type > 2 || args->type < 0)
        {
            fprintf(stdout, "type error input %d\n", args->type);
            print_help();
            exit(-1);
        }
    }
    if (argc >= 6)
    {
        args->sys_type = atoi(argv[5]);
        if (args->sys_type > 1 || args->sys_type < 0)
        {
            fprintf(stdout, "sys_type error input %d\n", args->sys_type);
            print_help();
            exit(-1);
        }
    }
    if (argc >= 7)
    {
        args->num_b = atoi(argv[6]);
        if (args->num_b < 0)
        {
            fprintf(stdout, "number of b error input %d\n", args->num_b);
            print_help();
            exit(-1);
        }
    }
    if (argc < 3 || argc > 7)
    {
        print_help();
        exit(-1);
    }
    if (argc < 6)
        read_mtx_header(args->filename_matrix.c_str(), &args->sys_type);

    fprintf(stdout, "matrix name :      %s\nvectorb name :     %s\nread function :    base-%d\ntype :             %d\nsys_type :         %d\nnum_b :          %d\n",
            args->filename_matrix.c_str(), args->filename_b.c_str(), args->read_matrix_base, args->type, args->sys_type, args->num_b);
}

int main(int argc, char **argv)
{
    long long int m, n, nnzA;
    int isSymmetricA;
    long long int *row_ptr = NULL;  // the csr row pointer array of matrix A
    long long int *col_idx = NULL;  // the csr column index array of matrix A
    double *val = NULL;             // the csr value array of matrix A (real number)
    double *val_im = NULL;          // the csr value array of matrix A (imaginary number)
    double *x = NULL, *x_im = NULL; // solution vector x, (x: real number, x_im: imaginary number)
    double *b = NULL, *b_im = NULL; // right-hand side vector b, (b: real number, b_im: imaginary number)
    double tt, time;

    SolverArgument args;

    /* ========================================== */
    // Step 0: Read command line argument
    /* ========================================== */

    //! If the user does not specify parameter sys_type,
    //! it will automatically determine whether it is a real system or a complex system according to the mtx file

    parse_args(argc, argv, &args);

    /* ========================================== */
    // Step 1: Load matrix and rhs from mtx files
    /* ========================================== */
    if (args.sys_type == 0) // real system
    {
        // load matrix
        int err = mmio_allinone_ilp64(&m, &n, &nnzA, &isSymmetricA, &args.read_matrix_base, &row_ptr, &col_idx, &val, args.filename_matrix.c_str());
        if (err == -1)
        {
            printf("Cannot open matrix file %s.\n", args.filename_matrix.c_str());
            return -1;
        }
        if (m != n)
        {
            fprintf(stdout, "Invalid matrix size.\n");
            return 0;
        }

        x = new double[n * args.num_b];
        b = new double[n * args.num_b];

        // load right-hand side vector b
        if (args.num_b == 1)
        {
            err = load_vector(n, b, args.filename_b.c_str());
            if (err == -1)
            {
                printf("Cannot open vector file %s.\n", args.filename_b.c_str());
                return -1;
            }
        }
        else
        {
            for (int i = 0; i < args.num_b; i++)
            {
                std::string filename = get_index_filename(args.filename_b, i).c_str();
                err = load_vector(n, &b[i * n], filename.c_str());
                if (err == -1)
                {
                    printf("Cannot open vector file %s.\n", filename.c_str());
                    return -1;
                }
            }
        }

        // initial vector x
        memset(x, 0.0, sizeof(double) * n * args.num_b);
    }
    else
    { // complex system
        mmio_allinone_complex_ilp64(&m, &n, &nnzA, &isSymmetricA, &args.read_matrix_base, &row_ptr, &col_idx, &val, &val_im, args.filename_matrix.c_str());
        if (m != n)
        {
            fprintf(stdout, "Invalid matrix size.\n");
            return 0;
        }

        x = new double[n * args.num_b];
        x_im = new double[n * args.num_b];
        b = new double[n * args.num_b];
        b_im = new double[n * args.num_b];

        // load right-hand side vector b
        if (args.num_b == 1)
            load_b_complex(n, b, b_im, args.filename_b.c_str());
        else
        {
            for (int i = 0; i < n; i++)
            {
                load_b_complex(n, &b[i * n], &b_im[i * n], get_index_filename(args.filename_b, i).c_str());
            }
        }

        // initial vector x
        memset(x, 0.0, sizeof(double) * n * args.num_b);
        memset(x_im, 0.0, sizeof(double) * n * args.num_b);
    }

    /* ========================================== */
    // Step 2: Solve the linear system
    /* ========================================== */
    TimerRec time_rec;

    for (int i = 0; i < args.test_frequency; i++)
    {
        if (args.sys_type == 0)
        {
            struct DirectSolver_ILP64 mysolver;
            mysolver.num_b = args.num_b;
            tt = GetCurrentTime();
            direct_preprocess_ilp64(&mysolver, n, row_ptr, col_idx);
            time_rec.preprocess_time += GetCurrentTime() - tt;

            tt = GetCurrentTime();
            direct_analyze_ilp64(&mysolver, n, val);
            time_rec.analyze_time += GetCurrentTime() - tt;

            tt = GetCurrentTime();
            direct_solve_ilp64(&mysolver, n, x, b);
            time_rec.solve_time += GetCurrentTime() - tt;
            direct_release_ilp64(&mysolver);
        }
        else
        {
            printf("ILP64 Complex Solver is not implemented\n");
            // struct DirectComplexSolver mysolver;

            // tt = GetCurrentTime();
            // direct_preprocess_complex(&mysolver, n, row_ptr, col_idx);
            // time_rec.preprocess_time += GetCurrentTime() - tt;

            // tt = GetCurrentTime();
            // direct_analyze_complex(&mysolver, n, val, val_im);
            // time_rec.analyze_time += GetCurrentTime() - tt;

            // tt = GetCurrentTime();
            // direct_solve_complex(&mysolver, n, x, x_im, b, b_im);
            // time_rec.solve_time += GetCurrentTime() - tt;
            // direct_release_complex(&mysolver);
        }
    }

    fprintf(stdout, "------------------------------------------\n");

    /* ========================================== */
    // Step 3: Check time, memory and correctness
    /* ========================================== */
    // check time
    if (args.type == 0)
    {
        time = (time_rec.preprocess_time + time_rec.analyze_time + time_rec.solve_time) / (double)args.test_frequency;
        fprintf(stdout, "CHECK end to end time :         %12.6lf ms\n", time);
    }
    else if (args.type == 1)
    {
        time = (time_rec.analyze_time + time_rec.solve_time) / (double)args.test_frequency;
        fprintf(stdout, "CHECK solver + solve time :     %12.6lf ms\n", time);
    }
    else if (args.type == 2)
    {
        time = time_rec.solve_time / (double)args.test_frequency;
        fprintf(stdout, "CHECK solve time :              %12.6lf ms\n", time);
    }

    // check the memory
    mem_usage();

    // store x to a file
    if (args.num_b == 1)
    {
        std::string answer_x = "answer_x.rhs";
        if (args.sys_type == 0) // real rhs
            store_x(n, x, answer_x.c_str());
        else // complex rhs
            store_x_complex(n, x, x_im, answer_x.c_str());
    }
    else
    {
        std::string answer_x = "answer_x";
        for (int i = 0; i < args.num_b; i++)
        {
            if (args.sys_type == 0) // real rhs
                store_x(n, &x[i * n], get_index_filename(answer_x, i).c_str());
            else // complex rhs
                store_x_complex(n, &x[i * n], &x_im[i * n], get_index_filename(answer_x, i).c_str());
        }
    }

    // check the correctness
    if (args.sys_type == 0) // real system
    {
        for (int i = 0; i < args.num_b; i++)
        {
            printf("Check for No. %d Solution:\n", i);
            // 1) using double precision
            check_correctness_ilp64(n, row_ptr, col_idx, val, &x[n * i], &b[n * i]);
            // 2) using long double precision
            check_correctness_ld_d2ld_ilp64(n, row_ptr, col_idx, val, &x[n * i], &b[n * i]);
        }
    }
    else
    { // complex system
        for (int i = 0; i < args.num_b; i++)
        {
            // printf("Check for No. %d Solution:\n", i);
            // // 1) using double precision
            // check_correctness_complex(n, row_ptr, col_idx, val, val_im, &x[n * i], &x_im[n * i], &b[n * i], &b_im[n * i]);
            // // 2) using long double precision
            // check_correctness_complex_ld_d2ld(n, row_ptr, col_idx, val, val_im, &x[n * i], &x_im[n * i], &b[n * i], &b_im[n * i]);
        }
        free(val_im);
        free(x_im);
        free(b_im);
    }

    free(row_ptr);
    free(col_idx);
    free(val);
    delete[] x;
    delete[] b;

    return 0;
}
