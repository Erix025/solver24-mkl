#!/bin/bash
#SBATCH --job-name=solver  # 作业名
#SBATCH --output=output/result_%j.txt  # 使用%j代表作业号，输出文件名将包含作业号
#SBATCH --error=output/error_%j.txt    # 错误输出文件也按作业号命名
#SBATCH --ntasks=1              # 运行任务的数量
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00         # 作业运行的时间（hh:mm:ss）
#SBATCH --partition=q_intel_share    # 作业所使用的分区
module load amd/Anaconda
module load amd/python
python ./analysis.py
