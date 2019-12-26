
# module load gcc/7.3.0 cuda/10.0.130
#!/bin/bash
#SBATCH --job-name=DeviceQuery
#SBATCH --nodes=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=10:00

nvidia-smi
make test_v1
