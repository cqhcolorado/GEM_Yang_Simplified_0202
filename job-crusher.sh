#!/bin/bash
#SBATCH -A FUS123 
#SBATCH -J coupling_test
#SBATCH -o %x-%j.out
#SBATCH -t 01:00:00
#SBATCH -p batch
### #SBATCH --reservation=hack3
#SBATCH -N 1

#!/bin/bash 

mkdir -p matrix
mkdir -p out
mkdir -p dump

export OMP_NUM_THREADS=7
srun -n 8 --ntasks-per-node=8 --gpus-per-task=1 -c 7 --gpu-bind=closest ./gem_main &> run.out 

#export OMP_NUM_THREADS=28
#srun -n 6 --ntasks-per-node=6 -c 7 -l  ./gem_main >run.out 2> run.err &

