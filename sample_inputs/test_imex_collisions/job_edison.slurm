#!/bin/bash -l
#SBATCH -q regular
#SBATCH -N 4
#SBATCH -t 01:00:00
#SBATCH -L SCRATCH

srun -n 96 ../../exec/cogent.Linux.64.CC.ftn.OPT.MPI.ex imex_collisionsnew.in >output.txt
