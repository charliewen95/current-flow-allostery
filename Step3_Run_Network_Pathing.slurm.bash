#!/bin/bash
#SBATCH -J cxPaths
#SBATCH --get-user-env
#SBATCH --partition=gpus
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=300:00:00

export PATH=/cm/shared/apps\:$PATH
which python


sys='n14k2'
rep="rep3 rep4"
nproc=12
echo "Starting pathing computations for $sys $rep using $nproc threads on `date`"
parallel -j $nproc "echo python  Compute_Pathing.py --inputPath cx26_GB_Network_Database/cx26_GB_Network.db -q '*, Seqid_1 + 226 * ( Chain_1 - 1 ) AS Resid_1, Seqid_2 + ( Chain_2 - 1 ) * 226 AS Resid_2 FROM Networks WHERE system = \"{1}\" and rep = \"{2}\" and Frame = 1' --sourceNodeNames {3} --targetNodeNames {4} --weightFunction 'reciprocal_1.' 'reciprocal_.61' --weightColumns 'Betweenness' 'TOTAL' --NodeColumns 'Resid_1' 'Resid_2' --groupingColumns 'system' 'rep' 'Frame' --stoppingCriteria 'convergence_.001' --outputDatabase Pathing.db -v -vl 2 ; python  Compute_Pathing.py --inputPath cx26_GB_Network_Database/cx26_GB_Network.db -q '*, Seqid_1 + 226 * ( Chain_1 - 1 ) AS Resid_1, Seqid_2 + ( Chain_2 - 1 ) * 226 AS Resid_2 FROM Networks WHERE system = \"{1}\" and rep = \"{2}\"' --sourceNodeNames {3} --targetNodeNames {4} --weightFunction 'reciprocal_1.' 'reciprocal_.61' --weightColumns 'Betweenness' 'TOTAL' --NodeColumns 'Resid_1' 'Resid_2' --groupingColumns 'system' 'rep' 'Frame' --stoppingCriteria 'convergence_.001' --outputDatabase Pathing.db -v -vl 2  > pathingLogFiles/Compute_Pathing.{1}.{2}.source_{3}.target_{4}.log" ::: $sys ::: $rep ::: `seq --separator " " 14 226 $((6*226))` ::: `seq --separator " " 47 226 $((6*226))`
echo "Done running pathing on `date`"
