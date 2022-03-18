#functions to facilitate loading path data from WISP logfiles
def load_wispLog_paths(logFilePath):
    pathList=[]
    with open(logFilePath,'r') as logFile:
        foundPathStart=False
        line=logFile.readline()
        while line and not foundPathStart:
            if 'Output identified paths' in line:
                foundPathStart=True
                #print line
            line=logFile.readline()
        if not foundPathStart:
            sys.stderr.write("Error: End of log file reached before path section was found")
        else:
            iPath=-1
            while '#' in line:
                if 'Path' in line:
                    if iPath>=0:
                        pathList.append(tempPath)
                        tempPath=[]
                        iPath=iPath+1
                    else:
                        tempPath=[]
                        iPath=iPath+1
                elif ('Length' in line):
                    tempPath.append(line.split()[2])
                elif ('Nodes:' in line):
                    tokens=np.array(line.split(),dtype='|S')
                    for node in tokens[2:len(tokens):2]:
                        tempPath.append(node)
                line=logFile.readline()
    return pathList

def simple_wispNode_converter(wispNode):
    return int(wispNode.split('_')[len(wispNode.split('_'))-1])

def wispPath_to_nodeIndPath(wispPath,
                            wispNode_to_nodeInd_function=simple_wispNode_converter):
    nodeIndPath=np.zeros(len(wispPath)-1)
    for iNode,wispNode in enumerate(wispPath[1:]):
        nodeIndPath[iNode]=wispNode_to_nodeInd_function(wispNode)
    return np.array(nodeIndPath,dtype=int)

def wispPaths_to_pathData(wispPaths,
                          wispNode_to_nodeInd_function=simple_wispNode_converter):
    pathData={
        'lengths':[],
        'paths':[],
        'count':len(wispPaths)
    }
    for iPath,wispPath in enumerate(wispPaths):
        pathData['lengths'].append(float(wispPath[0]))
        pathData['paths'].append(
            wispPath_to_nodeIndPath(wispPath,
                                    wispNode_to_nodeInd_function))
    return pathData

def load_wispLog_pathData(logFilePath,
                          wispNode_to_nodeInd_function=simple_wispNode_converter):
    return wispPaths_to_pathData(
        load_wispLog_paths(logFilePath),
        wispNode_to_nodeInd_function)

def netMatDict_to_nodeDataTable(netMatDict,
                                keySep='.',keyColNames=None,
                                indexCols=None,nodeIndName='X',
                                nodeInds=None,lapFunc=(lambda x: 2)):
    nodeTables=[]
    if (keyColNames is None):
        keyNames=['Key_%g'%iPart for iPart,part in \
                  enumerate(list(netMatDict.keys())[0].split(keySep))]
    else:
        keyNames=keyColNames
    if indexCols is None:
        indCols=keyNames
    else:
        indCols=indexCols
        for matKey in list(netMatDict.keys()):
            pbar.set_description_str(matKey)
            keyParts=matKey.split(keySep)
            tempMatLap=matLap(netMatDict[matKey])
            nodeVals=np.array(np.matrix(tempMatLap).diagonal()).flatten()
            norm=np.array(
                [lapFunc(
                    np.concatenate([tempMatLap[iRow,0:(iRow-1)],
                                    tempMatLap[iRow,(iRow+1):]])
                 ) for iRow,row in enumerate(tempMatLap)])
            nodeVals=nodeVals/norm
            nNodes=len(nodeVals)
            if nodeInds is None:
                inds=np.arange(nNodes)
            tempTable=pd.DataFrame({
                nodeIndName:inds,
                'value':nodeVals
            })
            for iCol,keyCol in enumerate(keyNames):
                tempTable[keyCol]=[keyParts[iCol]]*nNodes
            nodeTables.append(copy.deepcopy(tempTable))
            pbar.update()
    nodeTable=pd.concat(nodeTables)
    nodeDataTable=nodeTable.pivot_table(
        index=indCols,columns=nodeIndName,
        values='value',aggfunc=np.mean,fill_value=0
    )
    nodeDataTable.columns=np.array(nodeDataTable.columns)
    nodeDataTable=nodeDataTable.reset_index()
    return nodeDataTable

def netMatDict_to_edgeDataTable(netMatDict,
                                keySep='.',keyColNames=None,
                                indexCols=None,edgeIndNames=['X','Y'],
                                nodeInds=None,sparse=True):
    edgeTables=[]
    if (keyColNames is None):
        keyNames=['Key_%g'%iPart for iPart,part in \
                  enumerate(list(netMatDict.keys())[0].split(keySep))]
    else:
        keyNames=keyColNames
    if indexCols is None:
        indCols=keyNames
    else:
        indCols=indexCols
        for matKey in list(netMatDict.keys()):
            pbar.set_description_str('%s'%matKey)
            keyParts=matKey.split('.')
            tempMat=netMatDict[matKey]
            if sparse:
                edgeInds=np.nonzero(tempMat)
            else:
                pairs=np.array([[ii,jj] \
                       for ii in np.arange(tempMat.shape[0]) \
                       for jj in np.arange(tempMat.shape[1])])
                edgeInds=(pairs[:,0],pairs[:,1])
            edgeVals=np.array(tempMat)[edgeInds]
            nEdges=len(edgeVals)
            tempTable=pd.DataFrame({
                edgeIndNames[0]:edgeInds[0],
                edgeIndNames[1]:edgeInds[1],
                'value':edgeVals
            })
            for iCol,keyCol in enumerate(keyNames):
                tempTable[keyCol]=[keyParts[iCol]]*nEdges
            edgeTables.append(copy.deepcopy(tempTable))
            pbar.update()
    edgeTable=pd.concat(edgeTables)
    edgeDataTable=edgeTable.pivot_table(
        index=indCols,columns=edgeIndNames,
        values='value',aggfunc=np.mean,fill_value=0)
    edgeDataTable.columns=edgeDataTable.columns.map(
        lambda x: '_'.join(['%g'%xi for xi in x]))
    edgeDataTable=edgeDataTable.reset_index()
    return edgeDataTable


