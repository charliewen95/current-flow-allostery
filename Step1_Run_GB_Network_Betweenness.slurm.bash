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

parallel -j 12 "echo working on {3};python gb_network/Compute_GB_Betweenness_Frame.py -indir {1} -outdir {2} -i {3}.csv -o {3}.Betweenness -windmap -s `seq -s ' ' 14 226 1356` -t `seq -s ' ' 47 226 1356` -c 'Resid_1' 'Resid_2' -v -vl 2 --writeNodeVector -ft > {2}/{3}.log" ::: $inputDir ::: $outputDir ::: `ls $inputDir | sed "s/.csv//g"`
