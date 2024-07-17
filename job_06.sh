#!/bin/bash
#SBATCH --job-name=solver  # 作业名
#SBATCH --output=output/result_%j.txt  # 使用%j代表作业号，输出文件名将包含作业号
#SBATCH --error=output/error_%j.txt    # 错误输出文件也按作业号命名
#SBATCH --ntasks=1              # 运行任务的数量
#SBATCH --cpus-per-task=96
#SBATCH --time=01:00:00         # 作业运行的时间（hh:mm:ss）
#SBATCH --partition=q_intel_gpu_nvidia_h20_2    # 作业所使用的分区
#SBATCH --exclusive
#SBATCH --mem=1200G

module load intel/intel_compiler/2024u1
MATRIX=../../dataset/solverchallenge24_06/solverchallenge24_06_A.mtx
VECTOR=../../dataset/solverchallenge24_06/solverchallenge24_06_b.rhs
# MATRIX=matrix.mtx
# VECTOR=rhs.rhs
echo $MATRIX $VECTOR
export OMP_PROC_BIND=true
export OMP_NUM_THREADS=96
export MKL_NUM_THREADS=96
export MKL_PARDISO_OOC_MAX_CORE_SIZE=1000000
./build/direct_solver $MATRIX $VECTOR 1 0 1 1 ./ordering/ordering_6.rhs
# cp answer_x_*.rhs output/answer_$SLURM_JOB_ID
# ./iterative_solver $MATRIX $VECTOR 1
# OMP_NUM_THREADS=64 MKL_PARDISO_OOC_MAX_CORE_SIZE=200000 MKL_NUM_THREADS=64 ./hybrid_solver $MATRIX $VECTOR 1
