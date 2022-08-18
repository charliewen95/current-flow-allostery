try:
    import argparse
except ImportError:
    raise ImportError("require argparse")

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

#Others
import networkx as nx
import sqlite3
from sqlite3 import Error
import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


#self defined functions
if __name__ == "__main__":
    from functions import pathing_calc as pt_calc
else:
    from .functions import pathing_calc as pt_calc


########################################################

if __name__ == "__main__":
    parser=argparse.ArgumentParser(
        description="Loads the specified GB interaction network and calculates the corresponding flow betweenness network."+\
                    " This tool makes use of Yen's Algorithm from the networkx package to iteratively compute paths until "+\
                    " the desired stopping criteria is met."
    )
    
    parser.add_argument(
        '-uc','--useCSV',action='store',nargs='?',default=False,const=True,
        help='Turn this flag on to read and write data using CSV files instead of SQL databases'
    )
    
    dataLoadingParameters=parser.add_argument_group("Input / Data Loading Parameters")
    dataLoadingParameters.add_argument(
        '-i','--inputPath',default='./Network_Database_Directory/Network_Database',dest='inputPath',action='store',
        help='Path to the database (or csv file if --useCSV is set) containing the network data to be analyzed'
    )
    dataLoadingParameters.add_argument(
        '-q','--selectionQueryStrings',nargs='*',action='store',dest="selectionQueryStrings",
        help='List of query strings to select entries from the interaction data that specify the networks'+\
             '\nto be analyzed. If multiple strings are provided, the results are concatenated using pandas.concat. ' +\
             'if "--useCSV" is set, this will be fed to pandas.DataFrame.query, otherwise these will be used as ' +\
             'SQL queries on the database and the results appended row wise using pandas.concat.'
    )
    
    columnSetupParameters=parser.add_argument_group("Data Column Setup Parameters")
    columnSetupParameters.add_argument(
        '-gc','--groupingColumns',nargs='*',default=None,
        help='list of columns that will be used to define / group individual networks. '+\
             'Default is "None" which will assume that the entire loaded data set is a single network.'
    )
    columnSetupParameters.add_argument(
        '-cc','--computeResids',dest='computeResids',nargs='?',default=False,const=True,
        help='Use this option if you have nodes stored in terms of seqid and chain instead of simulation resid. '+\
             'be sure to also set "resPerChain" and "nChains" parameters (defaults are 226 and 6) and "seqCols" and "ChainCols"'
    )
    columnSetupParameters.add_argument(
        '-sqc','--seqCols',default=['Seqid_1','Seqid_2'],
        help='Name of the columns containing seqid (used when Resids need to be calculated)'
    )
    columnSetupParameters.add_argument(
        '-sqs','--seqStart',default=0,
        help='update this if sequence numbering starts at something other than 0 and resids are being calculated'
    )
    columnSetupParameters.add_argument(
        '-chs','--chainStart',default=0,
        help='update this if sequence numbering stats at something other than 0 and resids are being calculated'
    )
    columnSetupParameters.add_argument(
        '-rs','--residStart',default=0,
        help='Update this if resids are being calculated and numbering needs to start at something other than 0'
    )
    columnSetupParameters.add_argument(
        '-chc','--chainCols',default=['Chain_1','Chain_2'],
        help='Names of the columns containging chain numbers (used when Resids need to be calculated)'
    )
    columnSetupParameters.add_argument(
        '-rpc','--resPerChain',default=226,
        help='Number of residues per chain (used when Resids need to be calculated)'
    )
    columnSetupParameters.add_argument(
        '-nc','--nChains',default=6,
        help='Number of chains (Used when Resids need to be calculated)'
    )
    columnSetupParameters.add_argument(
        '-c','--NodeColumns',default=['Resid_1','Resid_2'],dest='nodeColumns',nargs='*',action='store',
        help='Names of the columns containing the names of the interacting nodes for each interaction entry'+\
             '\nexactly two arguments should be given. If not only the first two entries will get used.'
    )
    columnSetupParameters.add_argument(
        '-e','--weightColumns',default=['Betweenness'],nargs='*',action='store',dest="weightColumns",
        help='Name of the column used as edge weights. This column should not contain negative values.'+\
             ' As a sanity check, it will be fed through the absolute value function. See the --weightFunction '+\
             'argument for a list of other functions available (user defined weighting functions not available yet).'+\
             'if more than one weight column is specified, pathing will be repeated for each weight column given.'+\
             'If a function other than absolute value is needed, it can be specified using the "weightFunction" '+\
             'parameter. If there are multiple weight columns, a different weight function can be used for each '+\
             'weight column by providing a corresponding list of weight functions (see below). '+\
             'If you wish to use the same weight column with differrent weight functions, simply repeat its name '+\
             'in this list and supply an appropriate list of weight functions. If too few weight functions are '+\
             'supplied, the last supplied weight function will be repeated as needed.' 
    )
    columnSetupParameters.add_argument(
        '-w','--weightFunction',dest="weightFunction",nargs="*",default=['abs_1'],action='store',
        help='Function to apply to the weight column to convert it to weights. Defaults to "abs_1" (e.g. 1*abs(value)), '+\
             ' Other options are "reciprocal_SCALE" (SCALE/abs(value) e.g. "reciprocal_.6" = .6/abs(value)), '+\
             '"boltz_KT" (-KT*ln(abs(value)), e.g. "boltz_.61" would give -.61*ln(abs(value))), '+\
             '"rboltz_KT" (-KT*ln(1/abs(value)), '+\
             '"uinv_SCALE" (unit inverse, i.e. weight = 1 - abs(weightColumn)/SCALE. E.g. "uinv_.6" = 1-abs(value)/.6 ), '+\
             '"uharm_SCALE" (unit harmonic, i.e. weight = 1 - 1/(SCALE * abs(value)), '+\
             '"exp_LAMBDA" (corresponds to e**(-abs(value)/LAMBDA) e.g. "exp_.6" = e**(-abs(value)/.6) ), '+\
             '"rexp_LAMBDA" (e**(-1./(LAMBDA*abs(value))) e.g. "rexp_.6" = e**(-1/(.6*abs(value))), '+\
             '"uexp_LAMBDA" (1 - e**(-abs(value)/LAMBDA)), '+\
             '"hexp_LAMBDA" (1 - e**(-1/(LAMBDA*abs(value)))), '+\
             ' and "gauss_SIGMA" (corresponds to a gaussian with given standard deviation. E.g. '+\
             '"gauss_.25" = e**(-abs(value)**2 / (.25)**2)), '+\
             'and "rgauss_SIGMA" (e**(-1./(abs(value)*SIGMA)**2). '+\
             ' In most cases values of zero will get masked (corresponding edges are omitted).'+\
             ' User defined functions are not avaialable yet, but may be implemented in the future.'+\
             ' If a particular type is needed, SQL can be used to generate the needed column in most cases.'+\
             ' If multiple weight columns are given, then a list of weight functions may be provided to allow'+\
             ' different weight functions to be used for each weight column. If there are too few entries, then'+\
             ' the last weight function listed will be used for any remaining columns.'
    )
    columnSetupParameters.add_argument(
        '-s','--sourceNodeNames',nargs='+',dest='sourceNodeNames',action='store',
        help='string to be fed to the pandas DataFrame.query function to collect a list of source node names'+\
             '\nif multiple entries are given, each will be fed and the results aggregated into a list of'+\
             '\nof unique node names'
    )
    columnSetupParameters.add_argument(
        '-t','--targetNodeNames',nargs='+',dest='targetNodeNames',action='store',
        help='string to be fed to the pandas DataFrame.query function to collect a list of target node names'+\
             '\nif multiple entries are given, each will be fed and the results aggregated into a list of'+\
             '\nof unique node names. Note: if there is any overlap. I.e. sourceNodes and targetNodes litst'+\
             '\ncontain some of the same node(s) then a warning will be thrown and the common nodes will be put'+\
             '\ninto the source node list. The engine will be set to "python" to allow flexible selections'
    )
    
    outputParameters=parser.add_argument_group("Output parameters")
    outputParameters.add_argument(
        '-sc','--stoppingCriteria',dest="stoppingCriteria",default=["convergence_.0001"],nargs='*',action='store',
        help='stopping criteria, should be a string formatted as "TYPE_THRESHOLD" where "TYPE" is the '+\
             'stopping criteria type ("convergence", "relativeDilation", "totalDilation", or "count").' +\
             'And "THRESHOLD" a number describing the relevant threshold. '+\
             'For convergence, this is the relative change in combined node scores when a new path is added.'+\
             ' For dilation, this is either the increase in path length relative to the shortest path '+\
             'i.e. max_length=min_length*(1+valueParameter) or it is the total increase in path length '+\
             'i.e. max_length=min_length+valueParameter. For count, this is the total number of desired paths. '+\
             'The default is "convergence_.0001". If multiple entries are given, pathing will be run for each entry.'
    )
    outputParameters.add_argument(
        '-mp','-maxPaths',dest='maxPaths',default=10000,
        help='Maximum number of paths that will be calculated, functions as a safety to prevent excessive iterations if '+\
             'stopping criteria is progressing too slowly'
    )
    outputParameters.add_argument(
        '-o','--outputNameBase',default="Pathing",action='store',dest="outputNameBase",
        help='Base of the table names (or filepaths if "useCSV" is set) to write output to.'+\
             'For now, only paths are written. Node and edge usage may be adde later.'+\
             'The paths will be saved outFileNameBase.Paths'
    )
    outputParameters.add_argument(
        '-so','--outputDatabase',default=None,action='store',dest="outputDatabase",
        help='use this argument to specify the path to a different database that will be used for output. Has no effect if'+\
             ' "useCSV" is set'
    )
    outputParameters.add_argument(
        '-wto','--writeTimeout',default=30,
        help='Time to wait for database write lock to clear (no effect if "useCSV" is set). Defaults to 120 seconds '+\
             '(sqlite default is 5 seconds).'
    )
    outputParameters.add_argument(
        '-mwa','--maxWriteAttempts',default=4,
        help='Number of times to retry writting to database if a lock error occurs. '+\
             'I.e. total maximum wait time = maxWriteAttempts * writeTimeout. Defaults to 3 attempts'
    )
    outputParameters.add_argument(
        '-fscp','--failsafeCSVpath',default='./Pathing.Failsafe',
        help='Path used for writting a csv file in the event that writting fails due to repeated lock errors.'+\
             ' The process id of this thread and the local time (time.localtime()) will be added '+\
             'along with the ".csv" extension. E.g. Pathing.Failsafe.1234.2020_9_9_2_48_29_2_253.csv'
    )
    
    outputParameters.add_argument(
        '-dryrun','--dryrun',nargs='?',const=True,default=False,action='store',dest="dryrun",
        help='Dont run anything, jsut print out input argument namespace and end program'
    )
    outputParameters.add_argument(
        '-v','--verbose',nargs='?',const=True,default=False,action='store',dest="verbose",
        help='controls printing of progress / information to stdout during run'
    )
    outputParameters.add_argument(
        '-vl','--verboseLevel',default=0,action='store',dest="verboseLevel",
        help='when verbose flag is given, controls the amount of detail printed'
    )
    
    args=parser.parse_args()
 
    if args.verbose or args.dryrun:
        print('Input arguments:',args)
    if not args.dryrun:
        verbose=args.verbose
        verboseLevel=int(args.verboseLevel)
        outputNameBase=args.outputNameBase
        inputPath=args.inputPath
        
        timingsDict={
            "initialization":np.zeros(2),
            "loading":np.zeros(2),
            "setup":np.zeros(2),
            "pathFinding":np.zeros(2),
    		"save":np.zeros(2),
            "total":np.zeros(2)
        }
        timingsDict['total'][0]=time.time()
        timingsDict['initialization'][0]=time.time()
        def functionBuilder(fstring):
            fname=fstring.split('_')[0]
            fval=float(fstring.split('_')[1])
            fdict={
                "abs":lambda x: np.abs(x)/fval,
                "reciprocal":lambda x: fval/np.abs(x),
                "boltz":lambda x: -fval*np.log(np.abs(x)),
                "rboltz":lambda x: -fval*np.log(1./np.abs(x)),
                "uinv":lambda x: 1 - fval*np.abs(x),
                "uharm":lambda x: 1 - fval/np.abs(x),
                "exp":lambda x: np.exp(-np.abs(x)/fval),
                "rexp":lambda x: np.exp(-1./(fval*np.abs(x))),
                "uexp":lambda x: 1-np.exp(-np.abs(x)/fval),
                "hexp":lambda x: 1-np.exp(-1./(fval*np.abs(x))),
                "gauss":lambda x: np.exp(-(np.abs(x)/fval)**2),
                "rgauss":lambda x: np.exp(-1./(fval*np.abs(x))**2),
                "ugauss":lambda x: 1.-np.exp(-(np.abs(x)/fval)**2),
                "hgauss":lambda x: 1.-np.exp(-1./(fval*np.abs(x))**2)
            }
            return(fdict[fname])
           
        weightFuns=[functionBuilder(fstring) for fstring in args.weightFunction]
    
        timingsDict['initialization'][1]=time.time()
        timingsDict['loading'][0]=time.time()
    
        if verbose:
            print('loading data',end='\n' if args.verboseLevel==0 else ",")
       
        #Load data from disk
        if not args.useCSV:
            ### This code adapted from ###
            ### https://writeonly.wordpress.com/2009/07/16/simple-read-only-sqlalchemy-sessions/
            ### used to ensure input SQL query string will effectively not have write access
            def abort_ro(*args,**kwargs):
                ''' the terrible consequences for trying 
                    to flush to the db '''
                print("No writing allowed, tsk!  We're telling mom!")
                return 
            
            def db_setup(connstring='sqlite:///'+args.inputPath,
                        readOnly=True,echo=(verbose and (verboseLevel > 2))):
                engine = create_engine(connstring, echo=echo,connect_args={'timeout':float(args.writeTimeout)})
                Session = sessionmaker(
                    bind=engine, 
                    autoflush=not readOnly, 
                    autocommit=not readOnly
                )
                session = Session()
                
                if readOnly:
                    session.flush = abort_ro   # now it won't flush!
    
                return session, engine
            ### ### ###
            
            readSession, readEngine=db_setup()
            if not (args.outputDatabase is None):
                writeSession, writeEngine = db_setup(
                    connstring='sqlite:///'+args.outputDatabase,readOnly=False)
            else:
                writeSession,writeEngine=db_setup(readOnly=False)
            
            if verbose:
                print("Loading data using input SQL query")
                sys.stdout.flush()
                sys.stdout.flush()
                if verboseLevel>0:
                    t1=time.time()
            
            networkTables=[]
            for querySQL in args.selectionQueryStrings:
                query = readSession.query(querySQL)
                networkTables.append(pd.read_sql(query.statement,readEngine))
            networkData=pd.concat(networkTables)
            networkTables=[]
            gc.collect()
            if args.computeResids:
                networkData[args.nodeColumns[0]]=(networkData[args.seqCols[0]]-int(args.seqStart))+\
                    int(args.resPerChain)* (
                            (networkData[args.chainCols[0]]-int(args.chainStart)) % int(args.nChains)
                        )+int(args.residStart)
                networkData[args.nodeColumns[1]]=(networkData[args.seqCols[1]]-int(args.seqStart))+\
                    int(args.resPerChain)* (
                            (networkData[args.chainCols[1]]-int(args.chainStart)) % int(args.nChains)
                        )+int(args.residStart)
    
        else:
            networkTables=[]
            networkTable=pd.read_csv(args.inputPath)
            for queryString in args.selectionQueryStrings:
                networkTables.append(networkTable.query(queryString).copy())
            networkData=pd.concat(networkTables)
            networkTables=[]
            networkTable=None
            gc.collect()
    
            if args.computeResids:
                networkData[args.nodeColumns[0]]=(networkData[args.seqCols[0]]-int(args.seqStart))+\
                    int(args.resPerChain)* (
                            (networkData[args.chainCols[0]]-int(args.chainStart)) % int(args.nChains)
                        )+int(args.residStart)
                networkData[args.nodeColumns[1]]=(networkData[args.seqCols[1]]-int(args.seqStart))+\
                    int(args.resPerChain)* (
                            (networkData[args.chainCols[1]]-int(args.chainStart)) % int(args.nChains)
                        )+int(args.residStart)
    
        timingsDict['loading'][1]=time.time()
        timingsDict['setup'][0]=time.time()
        if verbose and (verboseLevel>1):
            print('Loaded Network Data:')
            print(networkData.head())
            if verboseLevel>2:
                print(networkData.describe())
            sys.stdout.flush()
    
        if verbose and (verboseLevel>1):
            print('source node names:',args.sourceNodeNames)
            print('target node names:',args.targetNodeNames)
        
        nodeColumn_1,nodeColumn_2=args.nodeColumns
        if verbose and (verboseLevel > 0):
            print('Building node name to index maps')
        nodeNames=np.unique(np.sort(np.concatenate([
                networkData[nodeColumn_1].unique(),
                networkData[nodeColumn_2].unique()])))
        nameToIndTable=pd.DataFrame({
            'NodeNames':np.array(nodeNames,dtype=str),
            'NodeInds':np.arange(len(nodeNames))
        })
        if verbose and (verboseLevel>3):
            print(nameToIndTable)
    
        if verbose and (verboseLevel>0):
            print('building source node list')
        sourceNodeNames=np.array(args.sourceNodeNames)
        sourceNodes=np.array([
            nameToIndTable.set_index('NodeNames')['NodeInds'].loc[sourceNodeName] \
            for sourceNodeName in sourceNodeNames
        ])
        if verbose and (verboseLevel>1):
            print('source nodes:')
            print(pd.DataFrame({'NodeNames':sourceNodes,'MatrixIndices':sourceNodes}))
        
        if verbose and (verboseLevel>0):
            print('building target node list')
            sys.stdout.flush()
        targetNodeNames=np.array(args.targetNodeNames)
        targetNodes=np.array([
            nameToIndTable.set_index('NodeNames')['NodeInds'].loc[targetNodeName] \
            for targetNodeName in targetNodeNames
        ])
        if verbose and (verboseLevel>1):
            print('target nodes:')
            print(pd.DataFrame({'NodeNames':targetNodeNames,'MatrixIndices':targetNodes}))
    
        sourceGrid,targetGrid=np.meshgrid(sourceNodes,targetNodes)
        pathTables=[]
        timingsDict['pathFinding'][0]=time.time()
        if args.groupingColumns is None:
            netData=networkData
            for iWeight,weightColumn in enumerate(args.weightColumns):
                if iWeight < len(weightFuns):
                    weightFun=weightFuns[iWeight]
                    funName=args.weightFunction[iWeight]
                else:
                    weightFun=weightFuns[-1]
                    funName=args.weightFunction[-1]
    
                if verbose and (verboseLevel > 0):
                    print(weightColumn+'__'+funName,end="(Pair: ")
                netData['weight']=weightFun(netData[weightColumn].copy())
                netData=netData.dropna()
    
                netMat=np.array(sp.sparse.coo_matrix(
                    (netData['weight'],
                    (nameToIndTable.set_index('NodeNames')['NodeInds'].loc[netData[nodeColumn_1].map(str)],
                    nameToIndTable.set_index('NodeNames')['NodeInds'].loc[netData[nodeColumn_2].map(str)])),
                    shape=(len(nameToIndTable),len(nameToIndTable))
                ).todense())
                netGraph=nx.from_numpy_array(netMat)
                
                for stoppingCriteria in args.stoppingCriteria:
                    stopType,stopVal=stoppingCriteria.split('_')
                    if not (stopType in ['convergence','relativeDilation','totalDilation','count']):
                        print("unrecognized stopping criteria '{ctype}'".format(ctype=stopType) +\
                                "defaulting to 'convergence_.0001'")
                        stopType='convergence'
                        stopVal=.0001
    
                    for sourceInd,targetInd in zip(sourceGrid.flatten(),targetGrid.flatten()):
                        startTime=time.time()
                        if verbose and (verboseLevel > 0):
                            print(
                                nameToIndTable.set_index('NodeInds')['NodeNames'].loc[sourceInd]+'__'+\
                                nameToIndTable.set_index('NodeInds')['NodeNames'].loc[targetInd]
                            )
                        if verbose and (verboseLevel > 0):
                            print(stopType+stopVal,end=" ")
                        if stopType=='convergence':
                            stopVal=float(stopVal)
                            paths=pt_calc.converge_subopt_paths_betweenness(
                                inputNetwork=netGraph,source=sourceInd,target=targetInd,weight='weight',
                                maxPaths=int(args.maxPaths),tolerance=stopVal,giveAlphas=False,verbose=(verboseLevel>2)
                            )
                        elif stopType=='relativeDilation':
                            stopVal=float(stopVal)
                            pathGenerator=nx.shortest_simple_paths(
                                G=netGraph, source=sourceInd, target=targetInd, k=1, weight='weight'
                            )
                            shortestPath=next(pathGenerator)
                            paths=[shortestPath]
                            minPathLen=pt_calc.calculatePathLength(
                                pathGraph=netGraph,path=paths[0],weight='weight'
                            )
                            pathLen=minPathLen
                            if verbose and (verboseLevel > 2):
                                print('Min_Path_Length {pl}'.format(pl=pathLen))
                            while( 
                                (((pathLen-minPathLen)/minPathLen) < stopVal) and \
                                (len(paths)<int(args.maxPaths)) 
                            ):
                                paths.append(next(pathGenerator))
                                pathLen=pt_calc.calculatePathLength(
                                    pathGraph=netGraph,path=paths[-1],weight='weight'
                                )
                                if verbose and (verboseLevel > 2):
                                    print(
                                        '{:.2e}'.format((pathLen-minPathLen)/minPathLen),
                                        end=' '
                                    )
                                    if verboseLevel>4:
                                        print(paths[-1])
    
                        elif stopType=='totalDilation':
                            stopVal=float(stopVal)
                            pathGenerator=nx.shortest_simple_paths(
                                G=netGraph, source=sourceInd, target=targetInd, k=1, weight='weight'
                            )
                            shortestPath=next(pathGenerator)
                            paths=[shortestPath]
                            minPathLen=pt_calc.calculatePathLength(
                                pathGraph=netGraph,path=paths[0],weight='weight'
                            )
                            pathLen=minPathLen
                            if verbose and (verboseLevel > 2):
                                print('Min_Path_Length {pl}'.format(pl=pathLen))
                            while( 
                                ((pathLen-minPathLen) < stopVal) and \
                                (len(paths)<int(args.maxPaths)) 
                            ):
                                paths.append(next(pathGenerator))
                                pathLen=pt_calc.calculatePathLength(
                                    pathGraph=netGraph,path=paths[-1],weight='weight'
                                )
                                if verbose and (verboseLevel > 2):
                                    print(
                                        '{:.2e}'.format((pathLen-minPathLen)),
                                        end=' '
                                    )
                                    if verboseLevel>4:
                                        print(paths[-1])
                        else:
                            stopVal=int(stopVal)
                            paths=pt_calc.k_shortest_paths(
                                G=netGraph, source=sourceInd, target=targetInd, k=stopVal, weight='weight'
                            )
    
                        pathTable=pd.DataFrame({
                            "Path_Rank":np.concatenate([[iPath]*len(path) for iPath,path in enumerate(paths)]),
                            "Path_NodeRank":np.concatenate([np.arange(len(path)) for path in paths]),
                            "Path_Length":np.concatenate(
                                [
                                    [pt_calc.calculatePathLength(
                                        pathGraph=netGraph,path=path,weight='weight'
                                    )]*len(path) \
                                    for path in paths
                                ]
                            )
                        })
                        for colName,colVal in zip(args.groupingColumns,netName):
                            pathTable[colName]=colVal
                        pathTable["Weight_Name"]=weightName
                        pathTable["Source"]=nameToIndTable.set_index('NodeInds')['NodeNames'].loc[sourceInd]
                        pathTable["Target"]=nameToIndTable.set_index('NodeInds')['NodeNames'].loc[targetInd]
                        pathTable["StoppingCriteria"]=stoppingCriteria
                        pathTable=pathTable[
                            np.concatenate([
                                pathTable.columns[3:],
                                pathTable.columns[:3]
                            ])
                        ]
                        pathTables.append(pathTable.copy())
                        gc.collect()
    
                        stopTime=time.time()
                        if verbose and (verboseLevel>0):
                            print("({pft:.3f} seconds)".format(pft=stopTime-startTime))
                            sys.stdout.flush()
        else:
            networkGroups=networkData.groupby(args.groupingColumns)
            if verbose:
                if verboseLevel > 0:
                    print("Detected {nNets} networks.".format(nNets=len(networkGroups)))
                print("Running Pathing")
            for netGroup in networkGroups:
                netName,netData=netGroup
                if verbose and verboseLevel>0:
                    print(
                        ', '.join([str(colName)+'__'+str(colVal) for colName,colVal in zip(args.groupingColumns,list(netName))]),
                        end=":"    
                    )
    
                for iWeight,weightColumn in enumerate(args.weightColumns):
                    if iWeight < len(weightFuns):
                        weightFun=weightFuns[iWeight]
                        funName=args.weightFunction[iWeight]
                    else:
                        weightFun=weightFuns[-1]
                        funName=args.weightFunction[-1]
    
                    weightName=weightColumn+'__'+funName
                    if verbose and (verboseLevel > 0):
                        print(weightName,end="(Pair: ")
                    netData['weight']=weightFun(netData[weightColumn].copy())
                    netData=netData.dropna()
    
                    netMat=np.array(sp.sparse.coo_matrix(
                        (netData['weight'],
                        (nameToIndTable.set_index('NodeNames')['NodeInds'].loc[netData[nodeColumn_1].map(str)],
                        nameToIndTable.set_index('NodeNames')['NodeInds'].loc[netData[nodeColumn_2].map(str)])),
                        shape=(len(nameToIndTable),len(nameToIndTable))
                    ).todense())
    
                    netGraph=nx.from_numpy_array(netMat)
    
                    for stoppingCriteria in args.stoppingCriteria:
                        stopType,stopVal=stoppingCriteria.split('_')
                        if not (stopType in ['convergence','relativeDilation','totalDilation','count']):
                            print("unrecognized stopping criteria '{ctype}'".format(ctype=stopType) +\
                                    "defaulting to 'convergence_.0001'")
                            stopType='convergence'
                            stopVal=.0001
    
                        for sourceInd,targetInd in zip(sourceGrid.flatten(),targetGrid.flatten()):
                            startTime=time.time()
                            if verbose and (verboseLevel > 0):
                                print(
                                    nameToIndTable.set_index('NodeInds')['NodeNames'].loc[sourceInd]+'__'+\
                                    nameToIndTable.set_index('NodeInds')['NodeNames'].loc[targetInd],
                                    end=" "
                                )
                            
                            if verbose and (verboseLevel > 0):
                                print(stopType+str(stopVal),end=" ")
                            if stopType=='convergence':
                                stopVal=float(stopVal)
                                paths=pt_calc.converge_subopt_paths_betweenness(
                                    inputNetwork=netGraph,source=sourceInd,target=targetInd,weight='weight',
                                    maxPaths=int(args.maxPaths),tolerance=stopVal,giveAlphas=False,verbose=(verboseLevel>2)
                                )
                            elif stopType=='relativeDilation':
                                stopVal=float(stopVal)
                                pathGenerator=nx.shortest_simple_paths(
                                    G=netGraph, source=sourceInd, target=targetInd, k=1, weight='weight'
                                )
                                shortestPath=next(pathGenerator)
                                paths=[shortestPath]
                                minPathLen=pt_calc.calculatePathLength(
                                    pathGraph=netGraph,path=paths[0],weight='weight'
                                )
                                pathLen=minPathLen
                                if verbose and (verboseLevel > 2):
                                    print('Min_Path_Length {pl}'.format(pl=pathLen))
                                while( 
                                    (((pathLen-minPathLen)/minPathLen) < stopVal) and \
                                    (len(paths)<int(args.maxPaths)) 
                                ):
                                    paths.append(next(pathGenerator))
                                    pathLen=pt_calc.calculatePathLength(
                                        pathGraph=netGraph,path=paths[-1],weight='weight'
                                    )
                                    if verbose and (verboseLevel > 2):
                                        print(
                                            '{:.2e}'.format((pathLen-minPathLen)/minPathLen),
                                            end=' '
                                        )
                                        if verboseLevel>4:
                                            print(paths[-1])
    
                            elif stopType=='totalDilation':
                                stopVal=float(stopVal)
                                pathGenerator=nx.shortest_simple_paths(
                                    G=netGraph, source=sourceInd, target=targetInd, k=1, weight='weight'
                                )
                                shortestPath=next(pathGenerator)
                                paths=[shortestPath]
                                minPathLen=pt_calc.calculatePathLength(
                                    pathGraph=netGraph,path=paths[0],weight='weight'
                                )
                                pathLen=minPathLen
                                if verbose and (verboseLevel > 2):
                                    print('Min_Path_Length {pl}'.format(pl=pathLen))
                                while( 
                                    ((pathLen-minPathLen) < stopVal) and \
                                    (len(paths)<int(args.maxPaths)) 
                                ):
                                    paths.append(next(pathGenerator))
                                    pathLen=pt_calc.calculatePathLength(
                                        pathGraph=netGraph,path=paths[-1],weight='weight'
                                    )
                                    if verbose and (verboseLevel > 2):
                                        print(
                                            '{:.2e}'.format((pathLen-minPathLen)),
                                            end=' '
                                        )
                                        if verboseLevel>4:
                                            print(paths[-1])
                            else:
                                stopVal=int(stopVal)
                                paths=pt_calc.k_shortest_paths(
                                    G=netGraph, source=sourceInd, target=targetInd, k=stopVal, weight='weight'
                                )
    
                            pathTable=pd.DataFrame({
                                "Path_Inds":np.concatenate([[iPath]*len(path) for iPath,path in enumerate(paths)]),
                                "Path_Steps":np.concatenate([np.arange(len(path)) for path in paths]),
                                "Path_Nodes":nameToIndTable.set_index('NodeInds')['NodeNames'].loc[np.concatenate(paths)],
                                "Path_Lengths":np.concatenate(
                                    [
                                        [pt_calc.calculatePathLength(
                                            pathGraph=netGraph,path=path,weight='weight'
                                        )]*len(path) \
                                        for path in paths
                                    ]
                                )
                            })
                            for colName,colVal in zip(args.groupingColumns,netName):
                                pathTable[colName]=colVal
                            pathTable["Weight_Name"]=weightName
                            pathTable["Source"]=nameToIndTable.set_index('NodeInds')['NodeNames'].loc[sourceInd]
                            pathTable["Target"]=nameToIndTable.set_index('NodeInds')['NodeNames'].loc[targetInd]
                            pathTable["StoppingCriteria"]=stoppingCriteria
                            pathTable=pathTable[
                                np.concatenate([
                                    pathTable.columns[3:],
                                    pathTable.columns[:3]
                                ])
                            ]
                            pathTables.append(pathTable.copy())
                            gc.collect()
    
                            stopTime=time.time()
                            if verbose and (verboseLevel>0):
                                print("({pft:.3f} seconds)".format(pft=stopTime-startTime))
                                sys.stdout.flush()
        timingsDict['pathFinding'][1]=time.time()
        pathingTable=pd.concat(pathTables)
        if verbose and (verboseLevel > 2):
            print(pathingTable.head())
        if verbose and (verboseLevel>0):
            print("saving paths",end=" ")
        timingsDict['save'][0]=time.time()
        
        if args.useCSV:
            pathingTable.to_csv(
                outputNameBase+'.Paths.csv'
            )
        else:
            writeSuccess=False
            tryCount=0
            while (not writeSuccess) and (tryCount < int(args.maxWriteAttempts)):
                try:
                    pathingTable.to_sql(
                        args.outputNameBase+'.Paths',
                        con=writeEngine,if_exists='append'
                    )
                    writeSuccess=True
                except Exception as e:
                    tryCount+=1
                    print("Failed to write results to database due to an error, {nleft} trys remaining".format(
                       nleft=int(args.maxWriteAttempts)-tryCount
                    ))
                    print(e.__doc__)
                    #print(e.message)
                    sys.stderr.write(e.__doc__)
                    #sys.stderr.write(e.message)
                else:
                    writeSuccess=True
    
            if not writeSuccess:
                outfilepath=args.failsafeCSVpath+'.'+str(os.getpid())+'.'+\
                    '_'.join(list(map(str,list(time.localtime())[:-1])))+'.csv'  
                print("WARNING! Failed to write data to database due to repeated lock failures")
                print("Writting data to failsafe CSV file instead: "+outfilepath)
                sys.stderr.write("WARNING! Failed to write data to database due to repeated lock failures")
                sys.stderr.write("Writting data to failsafe CSV file instead: "+outfilepath)
                pathingTable.to_csv(
                    outfilepath
                )
        timingsDict['save'][1]=time.time()
        timingsDict['total'][1]=time.time()
        if verbose:
            print("DONE")
            if verboseLevel>0:
                print("Timings:")
                for keyName in timingsDict:
                    print(keyName+':',timingsDict[keyName][1]-timingsDict[keyName][0],'seconds')
##############################

def compute_pathing(selectionQueryStrings,sourceNodeNames,targetNodeNames,useCSV=False,inputPath='./output_2/GB_Network.db',groupingColumns=None,computeResids=False,seqCols=['Seqid_1','Seqid_2'],seqStart=0,chainStart=0,residStart=0,chainCols=['Chain_1','Chain_2'],resPerChain=None,nChains=None,nodeColumns=['Resid_1','Resid_2'],weightColumns=['Betweenness'],weightFunction=['Betweenness'],stoppingCriteria=["convergence_.0001"],maxPaths=10000,outputNameBase="Pathing",outputDatabase=None,writeTimeout=30,maxWriteAttempts=4,failsafeCSVpath="./Pathing.Failsafe",dryrun=False,verbose=False,verboseLevel=0):
    ##### Defining the variables
    """
    Loads the specified GB interaction network and calculates the corresponding flow betweenness network. This tool makes use of Yen's Algorithm from the networkx package to iteratively compute paths until the desired stopping criteria is met.
    
    Default
    -------
    selectionQueryStrings	default         REQUIRED	
    sourceNodeNames		default 	REQUIRED
    targetNodeNames		default 	REQUIRED
    useCSV			default 	False
    inputPath			default		'./output_2/GB_Network.db'
    groupingColumns		default		None
    computeResids		default		False
    seqCols			default		['Seqid_1','Seqid_2']
    seqStart			default		0
    chainStart			default 	0
    residStart			default 	0
    chainCols			default 	['Chain_1','Chain_2']
    resPerChain			default		None
    nChains			default 	None
    nodeColumns			default 	['Resid_1','Resid_2']
    weightColumns		default		['Betweenness']
    weightFunction		default 	['abs_1']
    stoppingCriteria		default		["convergence_.0001"]
    maxPaths			default		10000
    outputNameBase		default	 	"Pathing"
    outputDatabase		default		None
    writeTimeout		default		30
    maxWriteAttempts		default		4
    failsafeCSVpath		default		'./Pathing.Failsafe'
    dryrun                      default		False
    verbose                     default		False
    verboseLevel                default		0

    Example
    -------
    

    Other Notes
    -----------
    sourceNodeNames & targetNodeNames format is :
    ['#','#','#','#']
    keep in mind if following format is used it will produce "key error"
    [#,#,#,#]

    """
####
    if resPerChain == None:
        print('resPerChain MISSING: SETTING TO DEFAULT TEST CASE VALUE 226')
        resPerChain = 226
    if nChains == None:
        print('nChains MISSING: SETTING TO DEFAULT TEST CASE VALUE 6')
        nChains = 6
    if groupingColumns == None:
        print("groupingColumns MISSING: SETTING TO DEFAULT TEST CASE VALUE 'system'")
        groupingColumns = 'system'

####Start of Code
    if verbose or dryrun:
        print('Input arguments:',\
" ",\
selectionQueryStrings,\
" ",\
sourceNodeNames,\
" ",\
targetNodeNames,\
" ",\
useCSV,\
" ",\
inputPath,\
" ",\
groupingColumns,\
" ",\
computeResids,\
" ",\
seqCols,\
" ",\
seqStart,\
" ",\
chainStart,\
" ",\
residStart,\
" ",\
chainCols,\
" ",\
resPerChain,\
" ",\
nChains,\
" ",\
nodeColumns,\
" ",\
weightColumns,\
" ",\
weightFunction,\
" ",\
stoppingCriteria,\
" ",\
maxPaths,\
" ",\
outputNameBase,\
" ",\
outputDatabase,\
" ",\
writeTimeout,\
" ",\
maxWriteAttempts,\
" ",\
failsafeCSVpath,\
" ",\
dryrun,\
" ",\
verbose,\
" ",\
verboseLevel\
)
    if not dryrun:
        verbose=verbose
        verboseLevel=int(verboseLevel)
        outputNameBase=outputNameBase
        inputPath=inputPath
        
        timingsDict={
            "initialization":np.zeros(2),
            "loading":np.zeros(2),
            "setup":np.zeros(2),
            "pathFinding":np.zeros(2),
    		"save":np.zeros(2),
            "total":np.zeros(2)
        }
        timingsDict['total'][0]=time.time()
        timingsDict['initialization'][0]=time.time()
        def functionBuilder(fstring):
            fname=fstring.split('_')[0]
            fval=float(fstring.split('_')[1])
            fdict={
                "abs":lambda x: np.abs(x)/fval,
                "reciprocal":lambda x: fval/np.abs(x),
                "boltz":lambda x: -fval*np.log(np.abs(x)),
                "rboltz":lambda x: -fval*np.log(1./np.abs(x)),
                "uinv":lambda x: 1 - fval*np.abs(x),
                "uharm":lambda x: 1 - fval/np.abs(x),
                "exp":lambda x: np.exp(-np.abs(x)/fval),
                "rexp":lambda x: np.exp(-1./(fval*np.abs(x))),
                "uexp":lambda x: 1-np.exp(-np.abs(x)/fval),
                "hexp":lambda x: 1-np.exp(-1./(fval*np.abs(x))),
                "gauss":lambda x: np.exp(-(np.abs(x)/fval)**2),
                "rgauss":lambda x: np.exp(-1./(fval*np.abs(x))**2),
                "ugauss":lambda x: 1.-np.exp(-(np.abs(x)/fval)**2),
                "hgauss":lambda x: 1.-np.exp(-1./(fval*np.abs(x))**2)
            }
            return(fdict[fname])
           
        weightFuns=[functionBuilder(fstring) for fstring in weightFunction]
    
        timingsDict['initialization'][1]=time.time()
        timingsDict['loading'][0]=time.time()
    
        if verbose:
            print('loading data',end='\n' if verboseLevel==0 else ",")
       
        #Load data from disk
        if not useCSV:
            ### This code adapted from ###
            ### https://writeonly.wordpress.com/2009/07/16/simple-read-only-sqlalchemy-sessions/
            ### used to ensure input SQL query string will effectively not have write access
            def abort_ro(*args,**kwargs):
                ''' the terrible consequences for trying 
                    to flush to the db '''
                print("No writing allowed, tsk!  We're telling mom!")
                return 
            
            def db_setup(connstring='sqlite:///'+inputPath,
                        readOnly=True,echo=(verbose and (verboseLevel > 2))):
                engine = create_engine(connstring, echo=echo,connect_args={'timeout':float(writeTimeout)})
                Session = sessionmaker(
                    bind=engine, 
                    autoflush=not readOnly, 
                    autocommit=not readOnly
                )
                session = Session()
                
                if readOnly:
                    session.flush = abort_ro   # now it won't flush!
    
                return session, engine
            ### ### ###
            
            readSession, readEngine=db_setup()
            if not (outputDatabase is None):
                writeSession, writeEngine = db_setup(
                    connstring='sqlite:///'+outputDatabase,readOnly=False)
            else:
                writeSession,writeEngine=db_setup(readOnly=False)
            
            if verbose:
                print("Loading data using input SQL query")
                sys.stdout.flush()
                sys.stdout.flush()
                if verboseLevel>0:
                    t1=time.time()
            
            networkTables=[]
###################################################
#            for querySQL in selectionQueryStrings:
#                query = readSession.query(querySQL)
#                networkTables.append(pd.read_sql(query.statement,readEngine))
            query = readSession.query(selectionQueryStrings)
            networkTables.append(pd.read_sql(query.statement,readEngine))
###################################################
            networkData=pd.concat(networkTables)
            networkTables=[]
            gc.collect()
            if computeResids:
                networkData[nodeColumns[0]]=(networkData[seqCols[0]]-int(seqStart))+\
                    int(resPerChain)* (
                            (networkData[chainCols[0]]-int(chainStart)) % int(nChains)
                        )+int(residStart)
                networkData[nodeColumns[1]]=(networkData[seqCols[1]]-int(seqStart))+\
                    int(resPerChain)* (
                            (networkData[chainCols[1]]-int(chainStart)) % int(nChains)
                        )+int(residStart)
    
        else:
            networkTables=[]
            networkTable=pd.read_csv(inputPath)
            for queryString in selectionQueryStrings:
                networkTables.append(networkTable.query(queryString).copy())
            networkData=pd.concat(networkTables)
            networkTables=[]
            networkTable=None
            gc.collect()
    
            if computeResids:
                networkData[nodeColumns[0]]=(networkData[seqCols[0]]-int(seqStart))+\
                    int(resPerChain)* (
                            (networkData[chainCols[0]]-int(chainStart)) % int(nChains)
                        )+int(residStart)
                networkData[nodeColumns[1]]=(networkData[seqCols[1]]-int(seqStart))+\
                    int(resPerChain)* (
                            (networkData[chainCols[1]]-int(chainStart)) % int(nChains)
                        )+int(residStart)
    
        timingsDict['loading'][1]=time.time()
        timingsDict['setup'][0]=time.time()
        if verbose and (verboseLevel>1):
            print('Loaded Network Data:')
            print(networkData.head())
            if verboseLevel>2:
                print(networkData.describe())
            sys.stdout.flush()
    
        if verbose and (verboseLevel>1):
            print('source node names:',sourceNodeNames)
            print('target node names:',targetNodeNames)
        
        nodeColumn_1,nodeColumn_2=nodeColumns
        if verbose and (verboseLevel > 0):
            print('Building node name to index maps')
        nodeNames=np.unique(np.sort(np.concatenate([
                networkData[nodeColumn_1].unique(),
                networkData[nodeColumn_2].unique()])))
        nameToIndTable=pd.DataFrame({
            'NodeNames':np.array(nodeNames,dtype=str),
            'NodeInds':np.arange(len(nodeNames))
        })
        if verbose and (verboseLevel>3):
            print(nameToIndTable)
    
        if verbose and (verboseLevel>0):
            print('building source node list')
        sourceNodeNames=np.array(sourceNodeNames)
        sourceNodes=np.array([
            nameToIndTable.set_index('NodeNames')['NodeInds'].loc[sourceNodeName] \
            for sourceNodeName in sourceNodeNames
        ])
        if verbose and (verboseLevel>1):
            print('source nodes:')
            print(pd.DataFrame({'NodeNames':sourceNodes,'MatrixIndices':sourceNodes}))
        
        if verbose and (verboseLevel>0):
            print('building target node list')
            sys.stdout.flush()
        targetNodeNames=np.array(targetNodeNames)
        targetNodes=np.array([
            nameToIndTable.set_index('NodeNames')['NodeInds'].loc[targetNodeName] \
            for targetNodeName in targetNodeNames
        ])
        if verbose and (verboseLevel>1):
            print('target nodes:')
            print(pd.DataFrame({'NodeNames':targetNodeNames,'MatrixIndices':targetNodes}))
    
        sourceGrid,targetGrid=np.meshgrid(sourceNodes,targetNodes)
        pathTables=[]
        timingsDict['pathFinding'][0]=time.time()
        if groupingColumns is None:
            netData=networkData
            for iWeight,weightColumn in enumerate(weightColumns):
                if iWeight < len(weightFuns):
                    weightFun=weightFuns[iWeight]
                    funName=weightFunction[iWeight]
                else:
                    weightFun=weightFuns[-1]
                    funName=weightFunction[-1]
    
                if verbose and (verboseLevel > 0):
                    print(weightColumn+'__'+funName,end="(Pair: ")
                netData['weight']=weightFun(netData[weightColumn].copy())
                netData=netData.dropna()
    
                netMat=np.array(sp.sparse.coo_matrix(
                    (netData['weight'],
                    (nameToIndTable.set_index('NodeNames')['NodeInds'].loc[netData[nodeColumn_1].map(str)],
                    nameToIndTable.set_index('NodeNames')['NodeInds'].loc[netData[nodeColumn_2].map(str)])),
                    shape=(len(nameToIndTable),len(nameToIndTable))
                ).todense())
                netGraph=nx.from_numpy_array(netMat)
                
                for stoppingCriteria in stoppingCriteria:
                    stopType,stopVal=stoppingCriteria.split('_')
                    if not (stopType in ['convergence','relativeDilation','totalDilation','count']):
                        print("unrecognized stopping criteria '{ctype}'".format(ctype=stopType) +\
                                "defaulting to 'convergence_.0001'")
                        stopType='convergence'
                        stopVal=.0001
    
                    for sourceInd,targetInd in zip(sourceGrid.flatten(),targetGrid.flatten()):
                        startTime=time.time()
                        if verbose and (verboseLevel > 0):
                            print(
                                nameToIndTable.set_index('NodeInds')['NodeNames'].loc[sourceInd]+'__'+\
                                nameToIndTable.set_index('NodeInds')['NodeNames'].loc[targetInd]
                            )
                        if verbose and (verboseLevel > 0):
                            print(stopType+stopVal,end=" ")
                        if stopType=='convergence':
                            stopVal=float(stopVal)
                            paths=corr_utils.converge_subopt_paths_betweenness(
                                inputNetwork=netGraph,source=sourceInd,target=targetInd,weight='weight',
                                maxPaths=int(maxPaths),tolerance=stopVal,giveAlphas=False,verbose=(verboseLevel>2)
                            )
                        elif stopType=='relativeDilation':
                            stopVal=float(stopVal)
                            pathGenerator=nx.shortest_simple_paths(
                                G=netGraph, source=sourceInd, target=targetInd, k=1, weight='weight'
                            )
                            shortestPath=next(pathGenerator)
                            paths=[shortestPath]
                            minPathLen=corr_utils.calculatePathLength(
                                pathGraph=netGraph,path=paths[0],weight='weight'
                            )
                            pathLen=minPathLen
                            if verbose and (verboseLevel > 2):
                                print('Min_Path_Length {pl}'.format(pl=pathLen))
                            while( 
                                (((pathLen-minPathLen)/minPathLen) < stopVal) and \
                                (len(paths)<int(maxPaths)) 
                            ):
                                paths.append(next(pathGenerator))
                                pathLen=corr_utils.calculatePathLength(
                                    pathGraph=netGraph,path=paths[-1],weight='weight'
                                )
                                if verbose and (verboseLevel > 2):
                                    print(
                                        '{:.2e}'.format((pathLen-minPathLen)/minPathLen),
                                        end=' '
                                    )
                                    if verboseLevel>4:
                                        print(paths[-1])
    
                        elif stopType=='totalDilation':
                            stopVal=float(stopVal)
                            pathGenerator=nx.shortest_simple_paths(
                                G=netGraph, source=sourceInd, target=targetInd, k=1, weight='weight'
                            )
                            shortestPath=next(pathGenerator)
                            paths=[shortestPath]
                            minPathLen=corr_utils.calculatePathLength(
                                pathGraph=netGraph,path=paths[0],weight='weight'
                            )
                            pathLen=minPathLen
                            if verbose and (verboseLevel > 2):
                                print('Min_Path_Length {pl}'.format(pl=pathLen))
                            while( 
                                ((pathLen-minPathLen) < stopVal) and \
                                (len(paths)<int(maxPaths)) 
                            ):
                                paths.append(next(pathGenerator))
                                pathLen=corr_utils.calculatePathLength(
                                    pathGraph=netGraph,path=paths[-1],weight='weight'
                                )
                                if verbose and (verboseLevel > 2):
                                    print(
                                        '{:.2e}'.format((pathLen-minPathLen)),
                                        end=' '
                                    )
                                    if verboseLevel>4:
                                        print(paths[-1])
                        else:
                            stopVal=int(stopVal)
                            paths=corr_utils.k_shortest_paths(
                                G=netGraph, source=sourceInd, target=targetInd, k=stopVal, weight='weight'
                            )
    
                        pathTable=pd.DataFrame({
                            "Path_Rank":np.concatenate([[iPath]*len(path) for iPath,path in enumerate(paths)]),
                            "Path_NodeRank":np.concatenate([np.arange(len(path)) for path in paths]),
                            "Path_Length":np.concatenate(
                                [
                                    [corr_utils.calculatePathLength(
                                        pathGraph=netGraph,path=path,weight='weight'
                                    )]*len(path) \
                                    for path in paths
                                ]
                            )
                        })
                        for colName,colVal in zip(groupingColumns,netName):
                            pathTable[colName]=colVal
                        pathTable["Weight_Name"]=weightName
                        pathTable["Source"]=nameToIndTable.set_index('NodeInds')['NodeNames'].loc[sourceInd]
                        pathTable["Target"]=nameToIndTable.set_index('NodeInds')['NodeNames'].loc[targetInd]
                        pathTable["StoppingCriteria"]=stoppingCriteria
                        pathTable=pathTable[
                            np.concatenate([
                                pathTable.columns[3:],
                                pathTable.columns[:3]
                            ])
                        ]
                        pathTables.append(pathTable.copy())
                        gc.collect()
    
                        stopTime=time.time()
                        if verbose and (verboseLevel>0):
                            print("({pft:.3f} seconds)".format(pft=stopTime-startTime))
                            sys.stdout.flush()
        else:
            networkGroups=networkData.groupby(groupingColumns)
            if verbose:
                if verboseLevel > 0:
                    print("Detected {nNets} networks.".format(nNets=len(networkGroups)))
                print("Running Pathing")
            for netGroup in networkGroups:
                netName,netData=netGroup
                if verbose and verboseLevel>0:
                    print(
                        ', '.join([str(colName)+'__'+str(colVal) for colName,colVal in zip(groupingColumns,list(netName))]),
                        end=":"    
                    )
    
                for iWeight,weightColumn in enumerate(weightColumns):
                    if iWeight < len(weightFuns):
                        weightFun=weightFuns[iWeight]
                        funName=weightFunction[iWeight]
                    else:
                        weightFun=weightFuns[-1]
                        funName=weightFunction[-1]
    
                    weightName=weightColumn+'__'+funName
                    if verbose and (verboseLevel > 0):
                        print(weightName,end="(Pair: ")
                    netData['weight']=weightFun(netData[weightColumn].copy())
                    netData=netData.dropna()
    
                    netMat=np.array(sp.sparse.coo_matrix(
                        (netData['weight'],
                        (nameToIndTable.set_index('NodeNames')['NodeInds'].loc[netData[nodeColumn_1].map(str)],
                        nameToIndTable.set_index('NodeNames')['NodeInds'].loc[netData[nodeColumn_2].map(str)])),
                        shape=(len(nameToIndTable),len(nameToIndTable))
                    ).todense())
    
                    netGraph=nx.from_numpy_array(netMat)
    
                    for stoppingCriteria in stoppingCriteria:
                        stopType,stopVal=stoppingCriteria.split('_')
                        if not (stopType in ['convergence','relativeDilation','totalDilation','count']):
                            print("unrecognized stopping criteria '{ctype}'".format(ctype=stopType) +\
                                    "defaulting to 'convergence_.0001'")
                            stopType='convergence'
                            stopVal=.0001
    
                        for sourceInd,targetInd in zip(sourceGrid.flatten(),targetGrid.flatten()):
                            startTime=time.time()
                            if verbose and (verboseLevel > 0):
                                print(
                                    nameToIndTable.set_index('NodeInds')['NodeNames'].loc[sourceInd]+'__'+\
                                    nameToIndTable.set_index('NodeInds')['NodeNames'].loc[targetInd],
                                    end=" "
                                )
                            
                            if verbose and (verboseLevel > 0):
                                print(stopType+str(stopVal),end=" ")
                            if stopType=='convergence':
                                stopVal=float(stopVal)
                                paths=corr_utils.converge_subopt_paths_betweenness(
                                    inputNetwork=netGraph,source=sourceInd,target=targetInd,weight='weight',
                                    maxPaths=int(maxPaths),tolerance=stopVal,giveAlphas=False,verbose=(verboseLevel>2)
                                )
                            elif stopType=='relativeDilation':
                                stopVal=float(stopVal)
                                pathGenerator=nx.shortest_simple_paths(
                                    G=netGraph, source=sourceInd, target=targetInd, k=1, weight='weight'
                                )
                                shortestPath=next(pathGenerator)
                                paths=[shortestPath]
                                minPathLen=corr_utils.calculatePathLength(
                                    pathGraph=netGraph,path=paths[0],weight='weight'
                                )
                                pathLen=minPathLen
                                if verbose and (verboseLevel > 2):
                                    print('Min_Path_Length {pl}'.format(pl=pathLen))
                                while( 
                                    (((pathLen-minPathLen)/minPathLen) < stopVal) and \
                                    (len(paths)<int(maxPaths)) 
                                ):
                                    paths.append(next(pathGenerator))
                                    pathLen=corr_utils.calculatePathLength(
                                        pathGraph=netGraph,path=paths[-1],weight='weight'
                                    )
                                    if verbose and (verboseLevel > 2):
                                        print(
                                            '{:.2e}'.format((pathLen-minPathLen)/minPathLen),
                                            end=' '
                                        )
                                        if verboseLevel>4:
                                            print(paths[-1])
    
                            elif stopType=='totalDilation':
                                stopVal=float(stopVal)
                                pathGenerator=nx.shortest_simple_paths(
                                    G=netGraph, source=sourceInd, target=targetInd, k=1, weight='weight'
                                )
                                shortestPath=next(pathGenerator)
                                paths=[shortestPath]
                                minPathLen=corr_utils.calculatePathLength(
                                    pathGraph=netGraph,path=paths[0],weight='weight'
                                )
                                pathLen=minPathLen
                                if verbose and (verboseLevel > 2):
                                    print('Min_Path_Length {pl}'.format(pl=pathLen))
                                while( 
                                    ((pathLen-minPathLen) < stopVal) and \
                                    (len(paths)<int(maxPaths)) 
                                ):
                                    paths.append(next(pathGenerator))
                                    pathLen=corr_utils.calculatePathLength(
                                        pathGraph=netGraph,path=paths[-1],weight='weight'
                                    )
                                    if verbose and (verboseLevel > 2):
                                        print(
                                            '{:.2e}'.format((pathLen-minPathLen)),
                                            end=' '
                                        )
                                        if verboseLevel>4:
                                            print(paths[-1])
                            else:
                                stopVal=int(stopVal)
                                paths=corr_utils.k_shortest_paths(
                                    G=netGraph, source=sourceInd, target=targetInd, k=stopVal, weight='weight'
                                )
    
                            pathTable=pd.DataFrame({
                                "Path_Inds":np.concatenate([[iPath]*len(path) for iPath,path in enumerate(paths)]),
                                "Path_Steps":np.concatenate([np.arange(len(path)) for path in paths]),
                                "Path_Nodes":nameToIndTable.set_index('NodeInds')['NodeNames'].loc[np.concatenate(paths)],
                                "Path_Lengths":np.concatenate(
                                    [
                                        [corr_utils.calculatePathLength(
                                            pathGraph=netGraph,path=path,weight='weight'
                                        )]*len(path) \
                                        for path in paths
                                    ]
                                )
                            })
                            for colName,colVal in zip(groupingColumns,netName):
                                pathTable[colName]=colVal
                            pathTable["Weight_Name"]=weightName
                            pathTable["Source"]=nameToIndTable.set_index('NodeInds')['NodeNames'].loc[sourceInd]
                            pathTable["Target"]=nameToIndTable.set_index('NodeInds')['NodeNames'].loc[targetInd]
                            pathTable["StoppingCriteria"]=stoppingCriteria
                            pathTable=pathTable[
                                np.concatenate([
                                    pathTable.columns[3:],
                                    pathTable.columns[:3]
                                ])
                            ]
                            pathTables.append(pathTable.copy())
                            gc.collect()
    
                            stopTime=time.time()
                            if verbose and (verboseLevel>0):
                                print("({pft:.3f} seconds)".format(pft=stopTime-startTime))
                                sys.stdout.flush()
        timingsDict['pathFinding'][1]=time.time()
        pathingTable=pd.concat(pathTables)
        if verbose and (verboseLevel > 2):
            print(pathingTable.head())
        if verbose and (verboseLevel>0):
            print("saving paths",end=" ")
        timingsDict['save'][0]=time.time()
        
        if useCSV:
            pathingTable.to_csv(
                outputNameBase+'.Paths.csv'
            )
        else:
            writeSuccess=False
            tryCount=0
            while (not writeSuccess) and (tryCount < int(maxWriteAttempts)):
                try:
                    pathingTable.to_sql(
                        outputNameBase+'.Paths',
                        con=writeEngine,if_exists='append'
                    )
                    writeSuccess=True
                except Exception as e:
                    tryCount+=1
                    print("Failed to write results to database due to an error, {nleft} trys remaining".format(
                       nleft=int(maxWriteAttempts)-tryCount
                    ))
                    print(e.__doc__)
                    #print(e.message)
                    sys.stderr.write(e.__doc__)
                    #sys.stderr.write(e.message)
                else:
                    writeSuccess=True
    
            if not writeSuccess:
                outfilepath=failsafeCSVpath+'.'+str(os.getpid())+'.'+\
                    '_'.join(list(map(str,list(time.localtime())[:-1])))+'.csv'  
                print("WARNING! Failed to write data to database due to repeated lock failures")
                print("Writting data to failsafe CSV file instead: "+outfilepath)
                sys.stderr.write("WARNING! Failed to write data to database due to repeated lock failures")
                sys.stderr.write("Writting data to failsafe CSV file instead: "+outfilepath)
                pathingTable.to_csv(
                    outfilepath
                )
        timingsDict['save'][1]=time.time()
        timingsDict['total'][1]=time.time()
        if verbose:
            print("DONE")
            if verboseLevel>0:
                print("Timings:")
                for keyName in timingsDict:
                    print(keyName+':',timingsDict[keyName][1]-timingsDict[keyName][0],'seconds')

