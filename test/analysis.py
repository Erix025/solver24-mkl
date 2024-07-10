import numpy as np
import scipy.io as sio
import scipy.sparse as sp
import matplotlib.pyplot as plt

def is_symmetric(matrix):
    if matrix.shape[0] != matrix.shape[1]:
        return False
    
    return (matrix != matrix.T).nnz == 0
matrix_file = '../../dataset/solverchallenge24_10/solverchallenge24_10_A.mtx'
sparse_matrix = sio.mmread(matrix_file)

# sparse_matrix_csr = sp.csr_matrix(sparse_matrix)

# plt.spy(sparse_matrix_csr, markersize=1)
# plt.title("visualized")
# output_file = '10.png'
# plt.savefig(output_file, dpi=300, bbox_inches='tight')
# print("Symmetric: ", is_symmetric(sparse_matrix))

