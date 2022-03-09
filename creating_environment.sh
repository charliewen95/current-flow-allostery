#!/bin/bash
##### CREATING ENVIRONMENT
conda create -n flowNetwork python=3.6 <<-EOF
y
y
EOF

eval "$(conda shell.bash hook)"
conda activate flowNetwork

##### Installing dependencies in environment
conda install -c conda-forge numpy <<-EOF
y
EOF

conda install -c conda-forge pandas <<-EOF
y
EOF

conda install -c conda-forge ambertools=21 compilers <<-EOF
y
EOF

conda install -c conda-forge scipy <<-EOF
y
EOF

conda install -c conda-forge sqlalchemy <<-EOF
y
EOF

conda install -c conda-forge networkx <<-EOF
y
EOF

conda install -c conda-forge dataclasses <<-EOF
y
EOF

conda install -c conda-forge parallel <<-EOF
y
EOF


