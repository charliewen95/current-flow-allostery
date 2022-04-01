def compute_pathing(useCSV,inputPath,selectionQueryStrings,groupingColumns,computeResids,seqCols,seqStart,chainStart,residStart,chainCols,resPerChain,nChains,nodeColumns,weightColumns,weightFunction,sourceNodeNames,targetNodeNames,stoppingCriteria,maxPaths,outputNameBase,outputDatabase,writeTimeout,maxWriteAttempts,failsafeCSVpath,dryrun,verbose,verboseLevel):
    ##### Defining the variables
    """
    Loads the specified GB interaction network and calculates the corresponding flow betweenness network. This tool makes use of Yen's Algorithm from the networkx package to iteratively compute paths until the desired stopping criteria is met.
    
    useCSV			default 	False
    inputPath			default		'./output_3'
    selectionQueryStrings	default         REQUIRED	
    groupingColumns		default		None
    computeResids		default		False
    seqCols			default		['Seqid_1','Seqid_2']
    seqStart			default		0
    chainStart			default 	0
    residStart			default 	0
    chainCols			default 	['Chain_1','Chain_2']
    resPerChain			default		226
    nChains			default 	6
    nodeColumns			default 	['Resid_1','Resid_2']
    weightColumns		default		['Betweenness']
    weightFunction		default 	['abs_1']
    sourceNodeNames		default 	REQUIRED
    targetNodeNames		default 	REQUIRED
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
    """
    if useCSV == None:
        useCSV = False
    if inputPath == None:
        inputPath = './output_3'
    if selectionQueryStrings == None:
        print('Selection Query Strings NOT GIVEN. REQUIRED.')
    if computeResids: 
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
                            paths=corr_utils.converge_subopt_paths_betweenness(
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
                            minPathLen=corr_utils.calculatePathLength(
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
                                (len(paths)<int(args.maxPaths)) 
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
                                paths=corr_utils.converge_subopt_paths_betweenness(
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
                                minPathLen=corr_utils.calculatePathLength(
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
                                    (len(paths)<int(args.maxPaths)) 
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
