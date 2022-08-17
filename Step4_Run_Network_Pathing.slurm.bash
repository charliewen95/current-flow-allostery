#!/bin/bash
#SBATCH -J curr-flo-allo
#SBATCH --get-user-env
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=300:00:00

export PATH=/cm/shared/apps\:$PATH
which python


sys='n14k2'
rep="rep1 rep2"
nproc=4
echo "Starting pathing computations for $sys $rep using $nproc threads on `date`"
parallel -j $nproc "echo python  current_flow_allostery/compute_pathing.py --inputPath output_3/GB_Betweenness_Bootstrapped_KS.db -q '*, Seqid_1 + 226 * ( Chain_1 - 1 ) AS Resid_1, Seqid_2 + ( Chain_2 - 1 ) * 226 AS Resid_2 FROM Networks WHERE system = \"{1}\" and rep = \"{2}\" and Frame = 1' --sourceNodeNames {3} --targetNodeNames {4} --weightFunction 'reciprocal_1.' 'reciprocal_.61' --weightColumns 'Betweenness' 'TOTAL' --NodeColumns 'Resid_1' 'Resid_2' --groupingColumns 'system' 'rep' 'Frame' --stoppingCriteria 'convergence_.001' --outputDatabase output_4/Pathing.db -v -vl 2 ; python  current_flow_allostery/compute_pathing.py --inputPath output_3/GB_Betweenness_Bootstrapped_KS.db -q '*, Seqid_1 + 226 * ( Chain_1 - 1 ) AS Resid_1, Seqid_2 + ( Chain_2 - 1 ) * 226 AS Resid_2 FROM Networks WHERE system = \"{1}\" and rep = \"{2}\"' --sourceNodeNames {3} --targetNodeNames {4} --weightFunction 'reciprocal_1.' 'reciprocal_.61' --weightColumns 'Betweenness' 'TOTAL' --NodeColumns 'Resid_1' 'Resid_2' --groupingColumns 'system' 'rep' 'Frame' --stoppingCriteria 'convergence_.001' --outputDatabase output_4/Pathing.db -v -vl 2  > pathingLogFiles/Compute_Pathing.{1}.{2}.source_{3}.target_{4}.log" ::: $sys ::: $rep ::: `seq --separator " " 14 226 $((6*226))` ::: `seq --separator " " 47 226 $((6*226))`

echo "Done running pathing on `date`"

#parallel -j $nproc "python <<- EOF
#import current_flow_allostery
#current_flow_allostery.compute_pathing(\
#'test',\
#{3},\
#{4},\
#useCSV=False,\
#inputPath='./output_2/network_database_name.db',\
#groupingColumns=None,\
#computeResids=False,\
#seqCols=['Seqid_1','Seqid_2'],\
#seqStart=0,\
#chainStart=0,\
#residStart=0,\
#chainCols=['Chain_1','Chain_2'],\
#resPerChain=226,\
#nChains=6,\
#nodeColumns=['Resid_1','Resid_2'],\
#weightColumns=['Betweenness'],\
#weightFunction=['abs_1'],\
#stoppingCriteria=['convergence_.0001'],\
#maxPaths=10000,\
#outputNameBase='Pathing',\
#outputDatabase='output_4/Pathing.db',\
#writeTimeout=30,\
#maxWriteAttempts=4,\
#failsafeCSVpath='./output_4/Pathing.Failsafe',\
#dryrun=False,\
#verbose=False,\
#verboseLevel=2\
#)
#EOF" ::: $sys ::: $rep ::: `seq --separator " " 14 226 $((6*226))` ::: `seq --separator " " 47 226 $((6*226))`
#
#echo "Done running pathing on `date`"
