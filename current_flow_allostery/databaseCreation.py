#numerical / data packages
try:
    import numpy as np
    np.set_printoptions(threshold=10)
except ImportError:
    raise ImportError("require numpy")
try:
    import pandas as pd
except ImportError:
    raise ImportError("require pandas")
try:
    import scipy
    import scipy.sparse
    import scipy as sp
except ImportError:
    raise ImportError("require scipy")

#utilities
import os
import sys
import gc
import copy
import time

#others
import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import sqlite3
from sqlite3 import Error
import gc
import pytraj as pt

#self defined functions
from .functions import database_mod as db_m


def databaseCreation(output_1='output_1',output_2='output_2',network_database_name='GB_Network.db'):
    """
    This function is the main function to generate Network Database in SQL.
    
    Default
    -------
    output_1				'output_1'
    output_2				'output_2'
    network_database_name		'GB_Network.db'
     
    Example
    -------    
    current_flow_allostery.databaseCreation(\
    '{1}',\
    '{2}',\
    '{3}.db',\
    )   
 
    Other notes
    -----------
    
    """
    #####Default Settings of the Variables
    
    #output_1 = networkDataDir
    #output_2 = network_database_directory
    networkDataDir = output_1
    network_database_directory = output_2
    
    datafilenameKeywords=[
        'csv','System','Replica','Frame','Betweenness','EdgeBetweenness'
    ]
    
    edgeTableNames=[
        filename for filename in os.listdir(networkDataDir) \
        if np.all([
            keyword in filename for keyword in datafilenameKeywords
        ])
    ]
    
    tempFrame=pd.read_csv(networkDataDir+'/'+edgeTableNames[0])
    tempFrameColumnNames=tempFrame.columns
    tempFrameColumnTypes=tempFrame.dtypes.map(str)
    print('column names:',tempFrameColumnNames)
    print('column types:',tempFrameColumnTypes)
    tempFrame.head()
    #####
    conn=db_m.create_connection(network_database_directory+'/'+network_database_name)
    
    tempFrame.set_index([
        'system','subtype','rep','Frame'
    ]).to_sql("Networks",conn,if_exists='replace')
    
    '''conn.execute(
        """
        create table Networks as
        select * from tempFrame
        """
    )'''
    #####
    datafilenameKeywords=[
        'csv','System','Replica','Frame','Betweenness','EdgeBetweenness'
    ]
    
    edgeTableNames=[
        filename for filename in os.listdir(networkDataDir) \
        if np.all([
            keyword in filename for keyword in datafilenameKeywords
        ])
    ]
    print(edgeTableNames)
    
    #initialize sqlalchemy engine
    engine = sqlalchemy.create_engine(
        'sqlite:///'+network_database_directory+'/'+network_database_name,
        echo=False)
    
    edgeTableName = np.sort(edgeTableNames)
     
    for iTable in edgeTableName:
        tempFrame=pd.read_csv(networkDataDir+'/'+iTable)
        tempFrame.to_sql('Networks',con=engine,if_exists='append',index=False)
        #conn.commit()
        gc.collect()
    
    print('Done')
    #####
    conn=db_m.create_connection(network_database_directory+'/'+network_database_name)
    conn.row_factory=sqlite3.Row
    
    cur=conn.cursor()
    ##########################################################################
    cur.execute("""
        SELECT * FROM Networks WHERE ((system='n14y2') AND (rep='rep1') AND (Frame=1) AND 
            ((Seqid_1=14) OR (Seqid_2=14)))
        """)
    #########################################################################^
    rows=cur.fetchall()
    testFrame=pd.DataFrame(rows,columns=rows[0].keys())
    conn.close()
    testFrame.head()
    #####
    nodefilenameKeywords=[
        'csv','System','Replica','Frame','Betweenness','NodeBetweenness'
    ]
    
    nodeTableNames=[
        filename for filename in os.listdir(networkDataDir) \
        if np.all([
            keyword in filename for keyword in nodefilenameKeywords
        ])
    ]
    
    engine = sqlalchemy.create_engine(
        'sqlite:///'+network_database_directory+'/'+network_database_name,
        echo=False)
   
    nodeTableName = np.sort(nodeTableNames)
     
    for iTable in nodeTableName:
        system,rep,frame=list(map(
                lambda x: x.split('__')[-1],
                iTable.split('.')[1:4]
        ))
        tempFrame=pd.read_csv(networkDataDir+'/'+iTable)
        tempFrame['system']=system
        tempFrame['rep']=rep
        tempFrame['Frame']=int(frame)
        tempFrame.to_sql('Network_Node_Betweenness',con=engine,if_exists='append')
        #conn.commit()
        gc.collect()
    
    print('Done') 
    #####
    conn=db_m.create_connection(network_database_directory+'/'+network_database_name)
    conn.row_factory=sqlite3.Row
    #########################################################Add variable 
    nodeList=np.concatenate([
        [14+ii*226 for ii in np.arange(6)],
        [47+ii*226 for ii in np.arange(6)]
    ])
    #########################################################^
    nodeSel=' OR '.join([
        'NodeName=%g'%nodeName \
        for nodeName in nodeList
    ])
    cur=conn.cursor()
    cur.execute("""
        SELECT system, NodeName, AVG(Betweenness) 
        FROM Network_Node_Betweenness 
        WHERE NOT (
                {nodesel}
            )
        GROUP BY system, NodeName
        ORDER BY Betweenness DESC
        """.format(nodesel=nodeSel))
    rows=cur.fetchall()
    testFrame=pd.DataFrame(rows,columns=rows[0].keys())
    conn.close()
    testFrame.head(n=10)
    #####
    ############################################################################
    visStruc=pt.load('structure_files/visualization_structure.pdb',top='structure_files/visualization_structure.parm7')
    ###########################################################################^
    nRes=visStruc.topology.n_residues
    ############################################################################
    nChains=6
    ###########################################################################^
    resPerChain=nRes/nChains
    
    caInds=visStruc.topology.atom_indices('@CA')
    caInfoTable=pd.DataFrame({
        'AtomIndex':caInds+1,
        'Name':[visStruc.topology.atom(iAtom).name for iAtom in caInds],
        'Type':[visStruc.topology.atom(iAtom).type for iAtom in caInds],
        'Element':[visStruc.topology.atom(iAtom).element for iAtom in caInds],
        'Charge':[visStruc.topology.atom(iAtom).charge for iAtom in caInds],
        'Radius':[visStruc.topology.atom(iAtom).gb_radius for iAtom in caInds],
        'Resid':[int(visStruc.topology.atom(iAtom).resid+1) for iAtom in caInds],
        'Resname':[visStruc.topology.atom(iAtom).resname for iAtom in caInds],
        'Seqid':[int((visStruc.topology.atom(iAtom).resid % resPerChain)+1) for iAtom in caInds],
        'Chain':[int(visStruc.topology.atom(iAtom).chain+1) for iAtom in caInds],
        'X':visStruc.xyz[0,caInds,0],
        'Y':visStruc.xyz[0,caInds,1],
        'Z':visStruc.xyz[0,caInds,2]
    })
    ############################################################################
    caInfoTable['system']='wt2'
    ############################################################################
    caInfoTable.head()
    
    caInfoTable.to_sql('Alpha_Carbon_Structure_Data',con=engine,if_exists='append')

###OLD SELF CLEAN UP: PLEASE IGNORE
##################################################
#    db_m.create_new_db(output_2+'/'+network_database_name)
#    
#    datafilenameKeywords=[
#        'csv','System','Replica','Frame','Betweenness','EdgeBetweenness'
#    ]
#    
#    edgeTableNames=[
#        filename for filename in os.listdir(output_1) \
#        if np.all([
#            keyword in filename for keyword in datafilenameKeywords
#        ])
#    ]
#    
#    tempFrame=pd.read_csv(output_1+'/'+edgeTableNames[0])
#    tempFrameColumnNames=tempFrame.columns
#    tempFrameColumnTypes=tempFrame.dtypes.map(str)
#    tempFrame.head()
#
#    #only need to use this if you want to delete a whole table
#    
#    #if not (conn is None):
#    #    conn.close()
#    try:
#        conn=db_m.create_connection(output_2+'/'+network_database_name)
#        cur=conn.cursor()
#        cur.execute('DROP TABLE Network_Node_Betweenness')
#        conn.commit()
#        conn.close()
#        print('ERROR: Network_Node_Betweenness already exist')
#    except:
#        print('Network_Node_Betweenness does not exist')
#    else:
#        return
#    try:
#        conn=db_m.create_connection(output_2+'/'+network_database_name)
#        cur=conn.cursor()
#        cur.execute('DROP TABLE Alpha_Carbon_Structure_Data')
#        conn.commit()
#        conn.close()
#        print('ERROR: Alpha_Carbon_Structure_Data already exist')
#        sys.exit(1)
#    except:
#        print('Alpha_Carbon_Structure_Data does not exist')
#    else:
#        return
#    try:
#        conn=db_m.create_connection(output_2+'/'+network_database_name)
#        cur=conn.cursor()
#        cur.execute('DROP TABLE Networks')
#        conn.commit()
#        conn.close()
#        print('ERROR: Networks already exist')
#    except:
#        print('Networks does not exist')
#    else:
#        return
#    print('no table to delete') 
#    #create new network table from the loaded frame
#    conn=db_m.create_connection(output_2+'/'+network_database_name)
#    
#    tempFrame.set_index([
#        'system','subtype','rep','Frame'
#    ]).to_sql("Networks",conn,if_exists='replace')
#    
#    '''
#    conn.execute("""CREATE TABLE Networks as SELECT * FROM tempFrame""")
#    '''
#    
#    #create list of network csv data files
#    datafilenameKeywords=[
#        'csv','System','Replica','Frame','Betweenness','EdgeBetweenness'
#    ]
#    
#    edgeTableNames=[
#        filename for filename in os.listdir(output_1) \
#        if np.all([
#            keyword in filename for keyword in datafilenameKeywords
#        ])
#    ]
#    #initialize sqlalchemy engine
#    engine = sqlalchemy.create_engine(
#        'sqlite:///'+output_2+'/'+network_database_name,
#        echo=False)
#    edgeTableName=np.sort(edgeTableNames) 
#    for iTable in edgeTableName:
#        tempFrame=pd.read_csv(output_1+'/'+iTable)
#        try:
#            tempFrame.to_sql(Networks,con=engine,if_exists='append',index=False)
#        except OSError as exc:
#            if exc.errno == 36:
#                print('throw name too long')
#            else:
#                raise  # re-raise previously caught exception 
#        except IOError:
#            print('READ in ERROR NOTE: nothing is read')
#            return
#        except :
#            print(iTable)
#            #return
#        else:
#            print('read in success')
#        
#        #conn.commit()
#        gc.collect()
#     
#    print('edgeTableNames Done')
#    conn=db_m.create_connection(output_2+'/'+network_database_name)
#    conn.row_factory=sqlite3.Row
#    
#    cur=conn.cursor()
#    cur.execute("""
#        SELECT * FROM Networks WHERE ((system='wt2') AND (rep='rep1') AND (Frame=1) OR 
#        ((Seqid_1=14) OR (Seqid_2=14)))
#        """)
#    rows=cur.fetchall()
#    testFrame=pd.DataFrame(rows,columns=rows[0].keys())
#    conn.close()
#    testFrame.head()
#    
#    nodefilenameKeywords=[
#        'csv','System','Replica','Frame','Betweenness','NodeBetweenness'
#    ]
#    
#    nodeTableNames=[
#        filename for filename in os.listdir(output_1) \
#        if np.all([
#            keyword in filename for keyword in nodefilenameKeywords
#        ])
#    ]
#      
#    engine = sqlalchemy.create_engine(
#        'sqlite:///'+output_2+'/'+network_database_name,
#        echo=False)
#    nodeTableName=np.sort(nodeTableNames)
#    for iTable in nodeTableName:
#        system,rep,frame=list(map(
#                lambda x: x.split('__')[-1],
#                iTable.split('.')[1:4]
#        ))
#        tempFrame=pd.read_csv(output_1+'/'+iTable)
#        tempFrame['system']=system
#        tempFrame['rep']=rep
#        tempFrame['Frame']=int(frame)
#        tempFrame.to_sql('Network_Node_Betweenness',con=engine,if_exists='append')
#        #conn.commit()
#        gc.collect()
#    
#    print('nodeTableNames Done')      
#    
#    conn=db_m.create_connection(output_2+'/'+network_database_name)
#    conn.row_factory=sqlite3.Row
#    
#    nodeList=np.concatenate([
#        [14+ii*226 for ii in np.arange(6)],
#        [47+ii*226 for ii in np.arange(6)]
#    ])
#    nodeSel=' OR '.join([
#        'NodeName=%g'%nodeName \
#        for nodeName in nodeList
#    ])
#    cur=conn.cursor()
#    cur.execute("""
#        SELECT system, NodeName, AVG(Betweenness) 
#        FROM Network_Node_Betweenness 
#        WHERE NOT (
#                {nodesel}
#            )
#        GROUP BY system, NodeName
#        ORDER BY Betweenness DESC
#        """.format(nodesel=nodeSel))
#    rows=cur.fetchall()
#    testFrame=pd.DataFrame(rows,columns=rows[0].keys())
#    conn.close()
#    testFrame.head(n=10)
#    visStruc=pt.load('/data/charlie/wes_code_clean_up/current-flow-allostery/structure_files/visualization_structure.pdb',top='/data/charlie/wes_code_clean_up/current-flow-allostery/structure_files/visualization_structure.parm7')
#    
#    nRes=visStruc.topology.n_residues
#    nChains=6
#    resPerChain=nRes/nChains
#    
#    caInds=visStruc.topology.atom_indices('@CA')
#    caInfoTable=pd.DataFrame({
#        'AtomIndex':caInds+1,
#        'Name':[visStruc.topology.atom(iAtom).name for iAtom in caInds],
#        'Type':[visStruc.topology.atom(iAtom).type for iAtom in caInds],
#        'Element':[visStruc.topology.atom(iAtom).element for iAtom in caInds],
#        'Charge':[visStruc.topology.atom(iAtom).charge for iAtom in caInds],
#        'Radius':[visStruc.topology.atom(iAtom).gb_radius for iAtom in caInds],
#        'Resid':[int(visStruc.topology.atom(iAtom).resid+1) for iAtom in caInds],
#        'Resname':[visStruc.topology.atom(iAtom).resname for iAtom in caInds],
#        'Seqid':[int((visStruc.topology.atom(iAtom).resid % resPerChain)+1) for iAtom in caInds],
#        'Chain':[int(visStruc.topology.atom(iAtom).chain+1) for iAtom in caInds],
#        'X':visStruc.xyz[0,caInds,0],
#        'Y':visStruc.xyz[0,caInds,1],
#        'Z':visStruc.xyz[0,caInds,2]
#    })
#    caInfoTable['system']='wt2'
#    
#    caInfoTable.head()
#    
#    caInfoTable.to_sql('Alpha_Carbon_Structure_Data',con=engine,if_exists='append')
#    
#    network_database_name='testWrite1.db'
#    nodefilenameKeywords=[
#        'csv','System','Replica','Frame','Betweenness','NodeBetweenness'
#    ]
#    
#    nodeTableNames=[
#        filename for filename in os.listdir(output_1) \
#        if np.all([
#            keyword in filename for keyword in nodefilenameKeywords
#        ])
#    ]
#    
#    engine = sqlalchemy.create_engine(
#        'sqlite:///'+output_2+'/'+network_database_name,
#        echo=False)
#     
#    nTables=101
#    tables=[]
#    nodeTableName=np.sort(nodeTableNames)
#    for iTable in nodeTableName:
#        system,rep,frame=list(map(
#                lambda x: x.split('__')[-1],
#                iTable.split('.')[1:4]
#        ))
#        tempFrame=pd.read_csv(output_1+'/'+iTable)
#        tempFrame['system']=system
#        tempFrame['rep']=rep
#        tempFrame['Frame']=int(frame)
#        tables.append(tempFrame.copy())
#        if (len(tables)>=nTables):
#            tempFrame=pd.concat(tables)
#            tempFrame.to_sql('Network_Node_Betweenness',con=engine,if_exists='append')
#            tables=[]
#        #conn.commit()
#        gc.collect()
#    
#    print('Alpha_Carbon_Structure_Data Done')
#    
#    network_database_name='testWrite2.db'
#    nodefilenameKeywords=[
#        'csv','System','Replica','Frame','Betweenness','NodeBetweenness'
#    ]
#    
#    nodeTableNames=[
#        filename for filename in os.listdir(output_1) \
#        if np.all([
#            keyword in filename for keyword in nodefilenameKeywords
#        ])
#    ]
#    
#    engine = sqlalchemy.create_engine(
#        'sqlite:///'+output_2+'/'+network_database_name,
#        echo=False)
#    
#    nTables=101
#    tables=[]
#    nodeTableName=np.sort(nodeTableNames)
#    for iTable in nodeTableName:
#        system,rep,frame=list(map(
#                lambda x: x.split('__')[-1],
#                iTable.split('.')[1:4]
#        ))
#        tempFrame=pd.read_csv(output_1+'/'+iTable)
#        tempFrame['system']=system
#        tempFrame['rep']=rep
#        tempFrame['Frame']=int(frame)
#        tables.append(tempFrame.copy())
#        if (len(tables)>=nTables):
#            tempFrame=pd.concat(tables)
#            tempFrame.to_sql('Network_Node_Betweenness',con=engine,if_exists='append',
#                             method='multi')
#            tables=[]
#        #conn.commit()
#        gc.collect()
#    
#    print('testWrite2.db Done')      
#    
#    tempFrame[
#        tempFrame['NodeName'].isin(np.arange(5)).map(lambda x: not x)
#    ]
#    
#    
#    network_database_name='testWrite2.db'
#    engine = sqlalchemy.create_engine(
#        'sqlite:///'+output_2+'/'+network_database_name,
#        echo=False)
#    
#    Session = sessionmaker(
#        bind=engine
#    )
#    session = Session()
#    session.bulk_insert_mappings(
#        'MentorInformation', 
#        tempFrame.to_dict(orient="records"))
#    session.close()
#   
#    gbKSdir='../cx26_GB_KS_Testing/KS_Result_Tables/'
#    gbKSfiles=[filename for filename in os.listdir(gbKSdir) \
#               if 'Results_Summary.csv' in filename]
#    
#    engine=sqlalchemy.create_engine(
#        'sqlite:///'+network_database_name+'/GB_Betweenness_Bootstrapped_KS.db'
#    )
#    
#    gbKStable=pd.concat([
#        pd.read_csv(gbKSdir+filename) for filename in gbKSfiles
#    ])
#    gbKStable=gbKStable.rename(columns={
#        'RefDistName':'Reference_system',
#        'TestDistname':'Test_system',
#        'alpha':'Alpha',
#        'ref_Result':'Ref_Differs',
#        'test_Result':'Test_Differs',
#        'null_cut':'nullCut',
#        'ref_cut':'refCut',
#        'test_cut':'testCut'
#    })
#    gbKStable.head()
#    
#    gbKStable.to_sql(
#        "Edge_GB_KS_Results",
#        con=engine,if_exists='append'
#    )
    
    
