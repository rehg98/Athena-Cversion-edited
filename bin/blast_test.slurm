#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0:30:00
#SBATCH --mem-per-cpu=4G

#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=rgolant@princeton.edu

./athena -i ../tst/3D-mhd/athinput.blasttest

