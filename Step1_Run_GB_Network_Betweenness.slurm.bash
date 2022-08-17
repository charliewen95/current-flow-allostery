#!/bin/bash
#SBATCH -J curr-flow-allo
#SBATCH --get-user-env
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=300:00:00

export PATH=/cm/shared/apps\:$PATH
inputDir=Test_Data
outputDir=output_1

which python
python --version

#parallel -j 12 "echo working on {3};python current_flow_allostery/betweenness.py -indir {1} -outdir {2} -i {3}.csv -o {3}.Betweenness -windmap -s `seq -s ' ' 14 226 1356` -t `seq -s ' ' 47 226 1356` -c 'Resid_1' 'Resid_2' -v -vl 2 --writeNodeVector -ft > {2}/{3}.log" ::: $inputDir ::: $outputDir ::: `ls $inputDir | sed "s/.csv//g"`

parallel -j 4 "python <<-EOF
import current_flow_allostery as cfa
cfa.betweenness(\
'{1}',\
'{2}',\
'{3}.csv',\
'{3}.Betweenness',\
sourceNodeNames=['14','240','466','692','918','1144'],\
targetNodeNames=['47','273','499','725','951','1177'],\
writeFullTable=True,\
verboseLevel=2\
)
EOF" ::: $inputDir ::: 'output_1' ::: `ls $inputDir | sed "s/.csv//g"`

#python current_flow_allostery/betweenness.py -indir Test_Data -outdir output_1 -i EnergyData_Network.System__n14y2_acetyl.Replica__rep1.Frame__001.csv -o EnergyData_Network.System__n14y2_acetyl.Replica__rep1.Frame__001.Betweenness.csv -windmap -s `seq -s ' ' 14 226 1356` -t `seq -s ' ' 47 226 1356` -c 'Resid_1' 'Resid_2' -v -vl 2 --writeNodeVector -ft > output_1/EnergyData_Network.System__n14y2_acetyl.Replica__rep1.Frame__001.Betweenness.log

#python <<-EOF
#import current_flow_allostery
#current_flow_allostery.betweenness(\
#'Test_Data',\
#'output_2',\
#'EnergyData_Network.System__n14y2_acetyl.Replica__rep1.Frame__001.csv',\
#'EnergyData_Network.System__n14y2_acetyl.Replica__rep1.Frame__001.Betweenness.csv',\
#None,\
#['Resid_1','Resid_2'],\
#None,\
#['14','240','466','692','918','1144'],\
#['47','273','499','725','951','1177'],\
#True,\
#True,\
#True,\
#None,\
#True,\
#2\
#)
#EOF
