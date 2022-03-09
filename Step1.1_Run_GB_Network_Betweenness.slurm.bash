#!/bin/bash
#SBATCH -J cx26btw
#SBATCH --get-user-env
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=300:00:00

export PATH=/cm/shared/apps\:$PATH
inputDir=GB_Network_Data
outputDir=GB_BTW_Network_Data

which python
python --version

python gb_network/cx26_GB_DataBase_Testing.py
