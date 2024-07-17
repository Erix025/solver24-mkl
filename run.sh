MATRIX=../../dataset/solverchallenge24_08/solverchallenge24_08_A.mtx
VECTOR=../../dataset/solverchallenge24_08/solverchallenge24_08_b.rhs
# MATRIX=test/matrix.mtx
# VECTOR=test/rhs
export OMP_PROC_BIND=true
# OMP_NUM_THREADS=64 MKL_PARDISO_OOC_MAX_CORE_SIZE=200000 MKL_NUM_THREADS=16 ./build/direct_solver_ilp64 $MATRIX $VECTOR 1 0 0 
OMP_NUM_THREADS=96 MKL_PARDISO_OOC_MAX_CORE_SIZE=1000000 MKL_NUM_THREADS=96 ./build/direct_solver $MATRIX $VECTOR 1 1 1
# OMP_NUM_THREADS=64 MKL_PARDISO_OOC_MAX_CORE_SIZE=200000 MKL_NUM_THREADS=64 ./hybrid_solver $MATRIX $VECTOR 1
