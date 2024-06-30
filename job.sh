#!/bin/bash
#SBATCH -J "name"
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./MolekularDynmaik/build
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --mem=200MB
#SBATCH --get-user-env
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --mail-user=ge47jow@mytum.de
#SBATCH --time=00:30:00

./MolSim ../input/"name".xml 
