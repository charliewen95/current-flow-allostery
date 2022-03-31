Usage 
=====
By this step, you should have already installed the package.

This package is written in 2 ways:
   1. Specific for bash submission
   2. As a python function
        e.g.
        import current_flow_allostery
        current_flow_allostery.betweenness()

For further and more specific examples of how to use, check out the attached bash submission script. Within you will be able to identify 4 different types of submission:
   1. Parallel bash submission
   2. Parallel python function submission
   3. Single run bash submission
   4. Single run python function submission

Currently this package only support Step1, which is the betweenness calculation.

Next patch: fixing databaseTesting.py bugs --> saving the stored csv data into sqlite3 databases

For Testing
-----------
Data:
        In the main package, it should contain a folder named 'Test_Data' this folder contains sample data generated from GB, using connexin 26 and its mutant counterpart N14Y.

Simple test:
        1. python                               #starts python
        2. import current_flow_allostery        #imports the package
        3.                                      #calls the function to calculate betweenness
        current_flow_allostery.betweenness(\
        'Test_Data',\
        'output_1',\
        'EnergyData_Network.System__n14y2_acetyl.Replica__rep1.Frame__001.csv',\
        'EnergyData_Network.System__n14y2_acetyl.Replica__rep1.Frame__001.Betweenness.csv',\
        None,\
        ['Resid_1','Resid_2'],\
        None,\
        ['14','240','466','692','918','1144'],\
        ['47','273','499','725','951','1177'],\
        True,\
        True,\
        True,\
        None,\
        True,\
        2\
        )
        4. current_flow_allostery.databaseTesting('output_1','output_2','network_database_name.db') --> currently only specific for connexin

Things that may go wrong:
        - Check environment (with pip install, dependencies should be installed)
          if fails, consider using the attatched creating_environment.sh which creates a specific environment for this package
        - Wrong output folder name
        

