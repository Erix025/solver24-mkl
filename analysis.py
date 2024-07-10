import numpy as np
import scipy.io as sio
import scipy.sparse as sp
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.sparse.linalg import eigs

def is_structurally_symmetric(sparse_matrix):
    rows, cols = sparse_matrix.shape
    
    for i in tqdm(range(rows)):
        row_indices = sparse_matrix.indices[sparse_matrix.indptr[i]:sparse_matrix.indptr[i+1]]
        
        for idx in row_indices:
            if i >= cols:
                continue
            if not sparse_matrix[idx, i]:
                return False
            
    return True

def sparse_matrix_determinant(sparse_matrix):
    if sparse_matrix.shape[0] != sparse_matrix.shape[1]:
        raise ValueError()
    eigenvalues = eigs(sparse_matrix, k=min(6, sparse_matrix.shape[0]-1), which='SR')

    determinant = np.prod(eigenvalues)

    return determinant

matrix_file = '../../dataset/solverchallenge24_09/solverchallenge24_09_A.mtx'
sparse_matrix = sio.mmread(matrix_file).tocsr()

print(is_structurally_symmetric(sparse_matrix))

# sparse_matrix_csr = sp.csr_matrix(sparse_matrix)

# plt.spy(sparse_matrix_csr, markersize=1)
# plt.title("visualized")
# output_file = '10.png'
# plt.savefig(output_file, dpi=300, bbox_inches='tight')
# print("Symmetric: ", is_symmetric(sparse_matrix))
