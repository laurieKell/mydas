#!/bin/sh 

#SBATCH --time=00:20:00
#SBATCH --nodes=2
#SBATCH -A myproj_id
#SBATCH -p DevQ

module load intel/2019

mpirun -n 80 mpi-benchmarks/src/IMB-MPI1