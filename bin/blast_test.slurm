#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=0:30:00

#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=rgolant@princeton.edu


srun ./athena -i ../tst/3D-mhd/athinput.blasttest

