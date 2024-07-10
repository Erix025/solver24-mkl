MATRIX=../../dataset/solverchallenge24_03/solverchallenge24_03_A.mtx
VECTOR=../../dataset/solverchallenge24_03/solverchallenge24_03_b.rhs
# MATRIX=test/matrix.mtx
# VECTOR=test/rhs
export OMP_PROC_BIND=true
OMP_NUM_THREADS=48 MKL_PARDISO_OOC_MAX_CORE_SIZE=1000000 MKL_NUM_THREADS=48 ./build/direct_solver $MATRIX $VECTOR 1 0 0
# ./iterative_solver $MATRIX $VECTOR 1
# OMP_NUM_THREADS=64 MKL_PARDISO_OOC_MAX_CORE_SIZE=200000 MKL_NUM_THREADS=64 ./hybrid_solver $MATRIX $VECTOR 1