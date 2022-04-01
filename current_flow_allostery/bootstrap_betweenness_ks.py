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
    import scipy.stats
except ImportError:
    raise ImportError("require scipy (NOTE: check scipy submodules)")

#utilities
import sys
import gc
import time

#self defined functions
from .functions import database_mod as db_m

def bootstrap_betweenness(databasePath,outputDatabase,querySQL,systemGroupColumn,interactionGroupColumns,referenceSystems,valueColumn,alphas,testAllPairs,outputBase,outputToCSV,writeBootstrapDistributions,flushGroupCount,dryrun,verbose,verboseLevel):
    ##### Defining the variables
    """
    databasePath		INPUT MUST BE GIVEN
    ouputDatabase 		default databasePath
    querySQL      		default 'SELECT *; FROM Networks'
    systemGroupColumn		default 'system'
    interactionGroupColumns	default ['Seqid_1','Seqid_2','Chain_Delta']
    referenceSystems		default ['wt2']
    valueColumn			default 'Betweenness'
    alphas			default [.1]
    testAllPairs		default False
    outputBase			default 'Edge_Betweenness_KS'
    outputToCSV			default False
    writeBootstrapDistributions	default False
    flushGroupCount		default 1
    dryrun			default False
    verbose			default False
    verboseLevel		default 0
    """
    if databasePath == None:
        print('Path to the directory containing the interaction network file. (required)')
    if outputDatabase == None:
        outputDatabase = databasePath
    if querySQL == None:
        querySQL = 'SELECT *; FROM Networks'
    if systemGroupColumn == None:
        systemGroupColumn = 'system'
    if interactionGroupColumns == None:
        interactionGroupColumns = ['Seqid_1','Seqid_2','Chain_Delta']
    if referenceSystems == None:
        referenceSystems = ['wt2']
    if valueColumn == None:
        valueColumn = 'Betweenness'
    if alphas == None:
        alphas = [.1]
    if testAllPairs == None:
        testAllPairs = False
    if outputBase == None: 
        outputBase = 'Edge_Betweenness_KS'
    if outputToCSV == None:
        outputToCSV = False
    if writeBootstrapDistributions == None:
        writeBootstrapDistributions = False
    if flushGroupCount == None:
        flushGroupCount = 1
    if dryrun == None:
        dryrun = False
    if verbose == None:
        verbose = False
    if verboseLevel == None:
        verboseLevel = 0

    ##### Start of code 
    if verbose or dryrun:
        print('Input arguments:',databasePath,outputDatabase,querySQL,systemGroupColumn,interactionGroupColumns,referenceSystems,valueColumn,alphas,testAllPairs,outputBase,outputToCSV,writeBootstrapDistributions,flushGroupCount,dryrun,verbose,verboseLevel)
    if not dryrun:
        if verbose and verboseLevel>0:
            ti=time.time()
        if verbose:
            print("Starting up")
        ##### 
        connstring='sqlite:///'+databasePath
        
        echo=( verbose and (verboseLevel > 2))
         
        alphas=np.sort(np.unique(np.array(list(alphas))))
        readSession, readEngine= db_m.db_setup(connstring,True,echo)
        if not (outputDatabase is None):
            connstring='sqlite:///'+outputDatabase
            writeSession, writeEngine = db_m.db_setup(
                connstring,False,echo)
        else:
            writeSession,writeEngine= db_m.db_setup(connstring,False,echo)
        #####        
        if verbose:
            print("Loading data using input SQL query")
            sys.stdout.flush()
            if verboseLevel>0:
                t1=time.time()
        
        networkData = pd.read_sql(querySQL,readEngine)
        if len(networkData)==0:
            print('Input query yielded no data! Exiting')
        else:
            if verbose:
                if verboseLevel>0:
                    t2=time.time()
                    print('-Input data query time: %.2f seconds'%(t2-t1))
                    t1=time.time()
                    if verboseLevel>1:
                        print('-- Data Extracted From Query --')
                        print(networkData.head())
                        print(networkData.describe())
                        print('-- -- --')
                print('Splitting Reference and Test System Data')
                sys.stdout.flush()
            
            refData=networkData[
                networkData[systemGroupColumn].isin(referenceSystems)
            ]
            testData=networkData[
                networkData[systemGroupColumn].isin(referenceSystems).map(lambda x: not x)
            ]
            if verbose:
                if (verboseLevel>0):
                    t2=time.time()
                    print('-Splitting time: %.2f seconds'%(t2-t1))
                    if verboseLevel>1:
                        print('--- reference data ---')
                        print(refData.head())
                        print(refData.describe())
                        print('\n--- Test Data ---')
                        print(testData.head())
                        print(testData.describe())
                print("--- --- Running Bootstrapping --- ---")
                refSystemGroups=refData.groupby(systemGroupColumn)
                testSystemGroups=testData.groupby(systemGroupColumn)
                for refSystemGroup in refSystemGroups:
                    refSystemName,refSystemData=refSystemGroup
                    for testSystemGroup in testSystemGroups:
                        testSystemName,testSystemData=testSystemGroup
                        if verbose:
                            print('--- Testing',refSystemName,'vs',testSystemName,'---')
                        interactionGroups=refSystemData.groupby(interactionGroupColumns)
                        
                        resultTables=[]
                        bootstrapTables=[]
                        for iGroup,interactionGroup in enumerate(interactionGroups):
                            interactionName,refInteractionData=interactionGroup
                            interactionQuery=' and '.join([
                                "({colname} == {colval})".format(
                                    colname=cname,colval=cval
                                ) \
                                for cname,cval in zip(interactionGroupColumns,interactionName)
                            ])
                            if verbose:
                                print('-Testing: ',interactionQuery,end=": ")
                                if verboseLevel>0:
                                    t1=time.time()
                                sys.stdout.flush()
                            testInteractionData=testSystemData.query(interactionQuery)
                            refVals=np.array(refInteractionData[valueColumn])
                            testVals=np.array(testInteractionData[valueColumn])
                            jointVals=np.concatenate([refVals,testVals])
                            
                            nBootSamples=int(np.max([64,1./(np.min(alphas)/2.)**2]))
                            if verbose:
                                print('Bootstrapping Null',end=', ')
                                sys.stdout.flush()
                            nullBootData=np.zeros(nBootSamples)
                            for iBoot in np.arange(nBootSamples):
                                nullBootData[iBoot]=sp.stats.ks_2samp(
                                    np.random.choice(a=jointVals,size=len(jointVals),replace=True),
                                    jointVals
                                ).statistic
                            if verbose:
                                print('Ref',end=', ')
                                sys.stdout.flush()
                            refBootData=np.zeros(nBootSamples)
                            for iBoot in np.arange(nBootSamples):
                                refBootData[iBoot]=sp.stats.ks_2samp(
                                    np.random.choice(a=refVals,size=len(refVals),replace=True),
                                    jointVals
                                ).statistic
                            if verbose:
                                print('Test',end='; ')
                                sys.stdout.flush()
                            testBootData=np.zeros(nBootSamples)
                            for iBoot in np.arange(nBootSamples):
                                testBootData[iBoot]=sp.stats.ks_2samp(
                                    np.random.choice(a=testVals,size=len(testVals),replace=True),
                                    jointVals
                                ).statistic
                            
                            if writeBootstrapDistributions:
                                if verbose:
                                    print('Compiling Bootstrap Data Table',end="; ")
                                    sys.stdout.flush()
                                bootDataFrame=pd.DataFrame({
                                    'Null_KS':nullBootData,
                                    'Reference_KS':refBootData,
                                    'Test_KS':testBootData
                                })
                                bootDataFrame['Reference_'+systemGroupColumn]=refSystemName
                                bootDataFrame['Test_'+systemGroupColumn]=testSystemName
                                for gColName,gColVal in zip(interactionGroupColumns,interactionName):
                                    bootDataFrame[gColName]=gColVal
                                bootstrapTables.append(bootDataFrame.copy())
                                if (len(bootstrapTables)>=flushGroupCount) or \
                                   (iGroup == (len(interactionGroups)-1)):
                                    bootDataFrame=pd.concat(bootstrapTables)
                                    if verbose:
                                        print('Flushing Bootstrap Data',end='; ')
                                        sys.stdout.flush()
                                    if outputToCSV:
                                        bootDataFrame.to_csv(
                                            outputBase+'_Bootstrap_Data.csv')
                                    else:
                                        bootDataFrame.to_sql(
                                            outputBase+'_Bootstrap_Data',
                                            con=writeEngine,if_exists='append'
                                        )
                                    bootstrapTables=[]
                                    gc.collect()
                            
                            if verbose:
                                print('Compiling Results Table',end="; ")
                                sys.stdout.flush()
                            resultsFrame=pd.DataFrame({
                                'Alpha':alphas,
                                'nullCut':[
                                    np.quantile(nullBootData,q=1.-alpha/2.) \
                                    for alpha in alphas
                                ],
                                'refCut':[
                                    np.quantile(refBootData,q=alpha/2.) \
                                    for alpha in alphas
                                ],
                                'testCut':[
                                    np.quantile(testBootData,q=alpha/2.) \
                                    for alpha in alphas
                                ]
                            })
                            resultsFrame['Ref_Differs']=resultsFrame['refCut']>resultsFrame['nullCut']
                            resultsFrame['Test_Differs']=resultsFrame['testCut']>resultsFrame['nullCut']
                            resultsFrame['Reference_'+systemGroupColumn]=refSystemName
                            resultsFrame['Test_'+systemGroupColumn]=testSystemName
                            for gColName,gColVal in zip(interactionGroupColumns,interactionName):
                                resultsFrame[gColName]=gColVal
                            resultTables.append(resultsFrame.copy())
                            if (len(resultTables)>=flushGroupCount) or \
                               (iGroup == (len(interactionGroups)-1)):
                                if verbose:
                                    print('Flushing Results Tables',end=';')
                                    sys.stdout.flush()
                                resultsFrame=pd.concat(resultTables)
                                if outputToCSV:
                                    resultsFrame.to_csv(
                                        outputBase+'_Results.csv'
                                    )
                                else:
                                    resultsFrame.to_sql(
                                        outputBase+'_Results',
                                        con=writeEngine,if_exists='append'
                                    )
                                resultTables=[]
                                gc.collect()
                            
                            if verbose:
                                if verboseLevel>0:
                                    t2=time.time()
                                    print(' time=%.2f s'%(t2-t1),end="")
                                print("")
                                sys.stdout.flush()
                            gc.collect()
        if verbose and verboseLevel>0:
            tf=time.time()
            print('--- --- --- --- ---')
            print('Total Run Time: %.4f minutes'%((tf-ti)/60.))
