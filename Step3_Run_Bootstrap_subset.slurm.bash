#!/bin/bash
#SBATCH -J BootKS
#SBATCH --get-user-env
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=300:00:00

export PATH=/cm/shared/apps\:$PATH
which python

seqidList=`seq 193 1 226`

parallel -j 12 "echo python current_flow_allostery/bootstrap_betweenness_ks.py -dbp output_3/GB_Network.db -qs '* FROM Networks WHERE (Seqid_1={1})' -fgc 25 --outputDatabase output_3/GB_Betweenness_Bootstrapped_KS.db-wbd --alphas .05 .1 .15 -v -vl 2 ;date > output_3/Bootstrap_Betweenness_KS.Seqid_1__{1}.log; python current_flow_allostery/bootstrap_betweenness_ks.py -dbp output_3/GB_Network.db -qs '* FROM Networks WHERE (Seqid_1={1})' -fgc 25 --outputDatabase output_3/Betweenness_Bootstrapped_KS.db --alphas .05 .1 .15 -v -vl 2 >> output_3/Bootstrap_Betweenness_KS.Seqid_1__{1}.log;date >> output_3/Bootstrap_Betweenness_KS.Seqid_1__{1}.log" ::: $seqidList

echo "done"

#####Network
#parallel -j 4 "python <<- EOF
#import current_flow_allostery
#current_flow_allostery.bootstrap_betweenness(\
#'output_3/GB_Network.db',\
#'output_3/GB_Betweenness_Bootstrapped_KS.db-wbd',\
#querySQL='SELECT * FROM Networks WHERE (Seqid_1={1})',\
#alphas=[0.05, 0.1, 0.15],\
#writeBootstrapDistributions=True,\
#flushGroupCount=25,\
#verboseLevel=2\
#)
#EOF" ::: $seqidList

#####Betweenness_Bootstrapped_KS
#parallel -j 4 "python <<- EOF
#import current_flow_allostery
#current_flow_allostery.bootstrap_betweenness(\
#'output_3/GB_Network.db',\
#'output_3/GB_Betweenness_Bootstrapped_KS.db',\
#querySQL='SELECT * FROM Networks WHERE (Seqid_1={1})',\
#alphas=[0.05, 0.1, 0.15],\
#writeBootstrapDistributions=True,\
#flushGroupCount=25,\
#verboseLevel=2\
#)
#EOF" ::: $seqidList
#
#echo "done"
