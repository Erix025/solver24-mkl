#include <stdio.h>
#include <stdlib.h>

void generate_csr_matrix(int n, int** row_ptr, int** col_idx, double** values)
{
    // 示例稀疏矩阵的非零元素
    // 对于n = 10，我们可以生成一个对角线为主的稀疏矩阵
    int nnz = 3 * n - 2; // 对角线+次对角线

    *row_ptr = (int*)malloc((n + 1) * sizeof(int));
    *col_idx = (int*)malloc(nnz * sizeof(int));
    *values  = (double*)malloc(nnz * sizeof(double));

    int index = 0;
    for (int i = 0; i < n; i++) {
        (*row_ptr)[i] = index;

        if (i > 0) {
            (*col_idx)[index] = i - 1;
            (*values)[index]  = -1.0;
            index++;
        }

        (*col_idx)[index] = i;
        (*values)[index]  = 2.0;
        index++;

        if (i < n - 1) {
            (*col_idx)[index] = i + 1;
            (*values)[index]  = -1.0;
            index++;
        }
    }
    (*row_ptr)[n] = index;
}

void write_mtx_file(const char* filename, int n, int nnz, int* row_ptr, int* col_idx, double* values)
{
    FILE* f = fopen(filename, "w");
    if (!f) {
        printf("Cannot open file %s\n", filename);
        return;
    }

    fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%d %d %d\n", n, n, nnz);

    for (int i = 0; i < n; i++) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            fprintf(f, "%d %d %lf\n", i, col_idx[j], values[j]);
        }
    }

    fclose(f);
}

void write_rhs_file(const char* filename, int n)
{
    FILE* f = fopen(filename, "w");
    if (!f) {
        printf("Cannot open file %s\n", filename);
        return;
    }

    fprintf(f, "%%MatrixMarket matrix array real general\n");
    fprintf(f, "%d 1\n", n);

    for (int i = 0; i < n; i++) {
        fprintf(f, "%lf\n", 1.0);
    }

    fclose(f);
}

int main()
{
    int     n = 10;
    int *   row_ptr, *col_idx;
    double* values;

    generate_csr_matrix(n, &row_ptr, &col_idx, &values);

    // 计算非零元素个数
    int nnz = (n * 3) - 2;

    write_mtx_file("matrix.mtx", n, nnz, row_ptr, col_idx, values);
    write_rhs_file("rhs.rhs", n);

    free(row_ptr);
    free(col_idx);
    free(values);

    return 0;
}