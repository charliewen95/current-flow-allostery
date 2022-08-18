#!/bin/bash
#SBATCH -J curr-flo-allo
#SBATCH --get-user-env
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=300:00:00

export PATH=/cm/shared/apps\:$PATH
which python


sys='wt2'
rep='rep1'
nproc=4

#echo "Starting pathing computations for $sys $rep using $nproc threads on `date`"
#parallel -j $nproc "echo python  current_flow_allostery/compute_pathing.py --inputPath output_2/GB_Network.db -q '*, Seqid_1 + 226 * ( Chain_1 - 1 ) AS Resid_1, Seqid_2 + ( Chain_2 - 1 ) * 226 AS Resid_2 FROM Networks WHERE system = \"{1}\" and rep = \"{2}\" and Frame = 1' --sourceNodeNames {3} --targetNodeNames {4} --weightFunction 'reciprocal_1.' 'reciprocal_.61' --weightColumns 'Betweenness' 'TOTAL' --NodeColumns 'Resid_1' 'Resid_2' --groupingColumns 'system' 'rep' 'Frame' --stoppingCriteria 'convergence_.001' --outputDatabase output_4/Pathing.db -v -vl 2 ; python  current_flow_allostery/compute_pathing.py --inputPath output_2/GB_Network.db -q '*, Seqid_1 + 226 * ( Chain_1 - 1 ) AS Resid_1, Seqid_2 + ( Chain_2 - 1 ) * 226 AS Resid_2 FROM Networks WHERE system = \"{1}\" and rep = \"{2}\"' --sourceNodeNames {3} --targetNodeNames {4} --weightFunction 'reciprocal_1.' 'reciprocal_.61' --weightColumns 'Betweenness' 'TOTAL' --NodeColumns 'Resid_1' 'Resid_2' --groupingColumns 'system' 'rep' 'Frame' --stoppingCriteria 'convergence_.001' --outputDatabase output_4/Pathing.db -v -vl 2  > pathingLogFiles/Compute_Pathing.{1}.{2}.source_{3}.target_{4}.log" ::: $sys ::: $rep ::: `seq --separator " " 14 226 $((6*226))` ::: `seq --separator " " 47 226 $((6*226))`
#
#echo "Done running pathing on `date`"

parallel -j $nproc "python <<- EOF
import current_flow_allostery
current_flow_allostery.compute_pathing(\
'*, Seqid_1 + 226 * ( Chain_1 - 1 ) AS Resid_1, Seqid_2 + ( Chain_2 - 1 ) * 226 AS Resid_2 FROM Networks WHERE system = \"{1}\" and rep = \"{2}\" and Frame = 1;',\
['14','240','466','692','918','1144'],\
['47','273','499','725','951','1177'],\
groupingColumns=['system', 'rep', 'Frame'],\
weightFunction=['reciprocal_1.', 'reciprocal_.61'],\
weightColumns=['Betweenness', 'TOTAL'],\
stoppingCriteria=['convergence_.001'],\
outputDatabase='./output_4/Pathing.db',\
verbose=True,\
verboseLevel=2\
)
EOF" ::: $sys ::: $rep

echo "Done running pathing on `date`"

