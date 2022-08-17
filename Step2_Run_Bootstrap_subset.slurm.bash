#!/bin/bash
#SBATCH -J BootKS
#SBATCH --get-user-env
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=300:00:00

export PATH=/cm/shared/apps\:$PATH
which python

seqidList=`seq 193 1 226`

parallel -j 12 "python <<- EOF
import current_flow_allostery
current_flow_allostery.bootstrap_betweenness(\
'output_2/network_database_name.db',\
'output_2/GB_Betweenness_Bootstrapped_KS.db',\
'SELECT * FROM Networks WHERE (Seqid_1={1})',\
None,\
None,\
None,\
None,\
[0.05, 0.1, 0.15],\
False,\
None,\
False,\
True,\
25,\
False,\
True,\
2\
)
EOF" ::: $seqidList


