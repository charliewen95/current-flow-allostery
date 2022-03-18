#functions for running and parsing subopt files generated via the VMD 'subopt' command

def validate_subopt_file(filepath,verbose=False):
    if not os.path.exists(filepath):
        if verbose:
            print(filepath+" does not exist.")
        return False
    if not os.path.isfile(filepath):
        if verbose:
            print(filepath+" is not a file.")
        return False
    foundPathStart=False
    foundPathCount=False
    with open(filepath,'r') as suboptFile:
        for line in suboptFile:
            if re.search('The final paths are',line):
                foundPathStart=True
            if re.search('Number of paths is',line):
                foundPathCount=True
            if foundPathStart and foundPathCount:
                break
    if foundPathStart and foundPathCount:
        if verbose:
            print(filepath+' is a valid subopt file')
        return True
    else:
        if verbose:
            print(filepath+" : ", end=' ')
            if not foundPathStart:
                print("is missing paths section, ", end=' ')
            if not foundPathCount:
                print("is missing path count")
        return False

def get_subopt_pathCount(filepath,verbose=False):
    if not validate_subopt_file(filepath,verbose=verbose):
        print("ERROR! "+filepath+" is not a valid subopt file.")
        return -1
    else:
        with open(filepath,'r') as suboptFile:
            for line in suboptFile:
                if re.search('Number of paths is',line):
                    tokens=str.split(line)
                    pathCount=tokens[len(tokens)-1]
                    if verbose:
                        print('path count = '+str(pathCount))
                    if not str.isdigit(pathCount):
                        break
                    return pathCount
        print("ERROR! Something went wrong, the file seemed valid but had an invalid path count line")
        return -1

def get_subopt_pathData(filepath,verbose=False):
    pathData=collections.OrderedDict()
    if not validate_subopt_file(filepath,verbose=verbose):
        print("ERROR! "+filepath+" is not a valid subopt file.")
        return pathData
    else:
        foundPathStart=False
        foundPathCount=False
        pathData['paths']=[]
        pathData['lengths']=[]
        with open(filepath,'r') as suboptFile:
            for line in suboptFile:
                if re.search('Number of paths is',line):
                    foundPathCount=True
                    tokens=str.split(line)
                    pathCount=tokens[len(tokens)-1]
                    if verbose:
                        print('path count = '+str(pathCount))
                    if not str.isdigit(pathCount):
                        foundPathCount=False
                        pathCount=-1
                    break
                if not foundPathStart:
                    if re.search('The final paths are',line):
                        foundPathStart=True
                else:
                    tokens=list(map(int,re.sub('[,)(]','',line).split()))
                    tempPath=tokens[0:(len(tokens)-1)]
                    pathData['paths'].append(tempPath)
                    pathData['lengths'].append(tokens[len(tokens)-1])
        if not foundPathStart:
            print("ERROR! "+filepath+" seemed valid but path section was apparently absent")
            return pathData
        if not foundPathCount:
            print("Warning! final path count line was missing or ill-formed!")
            pathData['count']=len(pathData['paths'])
        else:
            if len(pathData['paths']) != int(pathCount):
                print("Warning! subopt file lists number of paths as", end=' ')
                print(str(pathCount), end=' ')
                print("but", end=' ')
                print(len(pathData['paths']), end=' ')
                print("paths were found.")
                pathData['count']=len(pathData['paths'])
            else:
                pathData['count']=pathCount
        return pathData

def run_external_subopt(networkMatFilePath,outputFilePath,
                        sourceNode,targetNode,dilationValue,
                        externalSuboptCommand='subopt',
                        returnSuboptData=False,
                        verbose=False):
    if sourceNode > targetNode:
        sID=targetNode
        tID=sourceNode
    else:
        sID=sourceNode
        tID=targetNode
    suboptCommand=' '.join([
            externalSuboptCommand,
            networkMatFilePath,
            outputFilePath,
            str(dilationValue),str(sID),str(tID)])
    if verbose:
        print('running external subopt command: '+suboptCommand)
    os.system(suboptCommand)
    if returnSuboptData:
        return get_subopt_pathData(outputFilePath+".out",verbose=verbose)
    

def run_subopt_till_pathCount(networkMatFilePath,sourceNode,targetNode,
                              outputDir,outputBaseName='subopt',
                              minPathCount=0,percentDilationIncrement=10,
                              externalSuboptCommand='subopt',
                              returnSuboptData=False,onlyFinalRun=True,
                              verbose=False,verboseLevel=0):
    subVerbose=(verbose and (verboseLevel > 0))
    if verbose:
        print('running iterative dilation till '+str(minPathCount)+' paths are attained:')
        if not subVerbose:
            print('(%dilation,pathCount):', end=' ')
    percentDilation=0
    dilationValue=0
    outputFileName='.'.join([
            outputBaseName,
            '_'.join([str(sourceNode),str(targetNode)]),
            '_'.join(['dilation',str(percentDilation)])])
    outputFilePath='/'.join([outputDir,outputFileName])
    run_external_subopt(networkMatFilePath,outputFilePath,
                        sourceNode,targetNode,dilationValue,
                        externalSuboptCommand,returnSuboptData=False,
                        verbose=subVerbose)
    suboptDataFilePath='.'.join([outputFilePath,'out'])
    if not validate_subopt_file(suboptDataFilePath):
        print("ERROR! subopt failed to generate a valid output file. Aborting")
        if returnSuboptData:
            if onlyFinalRun:
                return collections.OrderedDict()
            else :
                return [collections.OrderedDict()]
        else:
            return
    pathCount=get_subopt_pathCount(suboptDataFilePath)
    if verbose:
        print("(%g,%g)"%(float(percentDilation),float(pathCount)), end=' ')
        if not subVerbose:
            print(",", end=' ')
        else:
            print("")
    tempData=get_subopt_pathData(suboptDataFilePath,subVerbose)
    optPathLength=np.min(list(map(float,tempData['lengths'])))
    if returnSuboptData:
        if onlyFinalRun:
            suboptData=tempData
        else:
            suboptData=[tempData]
    while float(pathCount) < float(minPathCount):
        percentDilation=percentDilation+percentDilationIncrement
        dilationFactor=percentDilation/100.0
        dilationValue=optPathLength*dilationFactor
        outputFileName='.'.join([
            outputBaseName,
            '_'.join([str(sourceNode),str(targetNode)]),
            '_'.join(['dilation',str(percentDilation)])])
        outputFilePath='/'.join([outputDir,outputFileName])
        run_external_subopt(networkMatFilePath,outputFilePath,
                            sourceNode,targetNode,dilationValue,
                            externalSuboptCommand,returnSuboptData=False,
                            verbose=subVerbose)
        suboptDataFilePath='.'.join([outputFilePath,'out'])
        if not validate_subopt_file(suboptDataFilePath):
            print("ERROR! subopt failed to generate a valid output file. Aborting")
            if returnSuboptData:
                return suboptoptData
            else:
                return
        pathCount=get_subopt_pathCount(suboptDataFilePath)
        if verbose:
            print("(%g,%g)"%(float(percentDilation),float(pathCount)), end=' ')
        if not subVerbose:
            print(",", end=' ')
        else:
            print("")
        tempData=get_subopt_pathData(suboptDataFilePath,subVerbose)   
        if returnSuboptData:
            if onlyFinalRun:
                suboptData=tempData
            else:
                suboptData.append(tempData)
    if verbose and not subVerbose:
        print("")
    if verbose:
        print("DONE!")
    if returnSuboptData:
        return suboptData
    
def get_subopt_dilations_data(suboptDir,basePattern,sep='.',
                              onlyMaxDilation=True,includeDilationValues=False,
                              includeFileNames=False,
                              verbose=False,verboseLevel=0):
    fileSearchPattern=str(sep).join([basePattern,'_'.join(['dilation','*']),'out'])
    if verbose:
        print('file search pattern = '+fileSearchPattern)
    searchPathPattern='/'.join([suboptDir,fileSearchPattern])
    if verbose:
        print('searchPathPattern = '+searchPathPattern)
    filePathList=glob.glob(searchPathPattern)
    fileNameList=[filepath.split('/')[-1] for filepath in filePathList]
    if verbose and (verboseLevel > 0):
        print('file name list: '+'\n'.join(fileNameList))
    if not onlyMaxDilation:
        suboptDataSets=[]
    else:
        suboptData=collections.OrderedDict()
        maxDilationValue=-1
    for iFile in np.arange(len(filePathList)):
        suboptFilePath=filePathList[iFile]
        suboptFileName=fileNameList[iFile]
        if validate_subopt_file(suboptFilePath):
            suboptData=get_subopt_pathData(suboptFilePath)
            nameTokens=suboptFileName.split(sep)
            dilationToken=[token for token in nameTokens if 'dilation' in token][-1]
            dilationValue=int(dilationToken.split("_")[-1])
            if includeDilationValues:
                suboptData['dilation']=dilationValue
            if includeFileNames:
                suboptData['filename']=suboptFileName
            if onlyMaxDilation:
                if dilationValue > maxDilationValue:
                    suboptDataSets=suboptData
                    maxDilationValue=dilationValue
            else:
                suboptDataSets.append(suboptData)
        else:
            print("Warning: "+suboptFileName+" was not a valid subopt file.")
    if onlyMaxDilation:
        return suboptData
    else:
        return suboptDataSets

def get_index_of_nth_maxRank_element(valList,n,verbose=False):
    #this assumes that valList is already sorted in ascending order!
    u,v=np.unique(valList,return_inverse=True)
    maxRanks=(np.cumsum(np.bincount(v,minlength=u.size))-1)[v]
    if verbose:
        print('value array:', end=' ')
        print(valList)
        print('rank array:', end=' ')
        print(maxRanks)
    if np.max(maxRanks) < n:
        print("get nth maxRank: Warning! there are not enough elements")
        return(len(valList)-1)
    else:
        return [m for m,i in enumerate(maxRanks) if i >= n][0]

def get_top_n_pathData_paths(pathData,n,verbose=False,verboseLevel=0):
    subVerbose=(verbose and (verboseLevel>0))
    outData=copy.deepcopy(pathData)
    pathSortingArray=np.argsort(outData['lengths'])
    outData['paths']=[outData['paths'][i] for i in pathSortingArray]
    outData['lengths']=[outData['lengths'][i] for i in pathSortingArray]
    maxPathIndex=get_index_of_nth_maxRank_element(outData['lengths'],n,
                                                  verbose=subVerbose)
    if verbose:
        print('max path index = '+str(maxPathIndex))
    outData['paths']=outData['paths'][0:(maxPathIndex+1)]
    outData['lengths']=outData['lengths'][0:(maxPathIndex+1)]
    outData['count']=len(outData['lengths'])
    return outData

def get_pathData_node_count_array(pathData,nNodes,inputIndexBase=0,verbose=False):
    outArray=np.zeros(nNodes)
    for path in pathData['paths']:
        outArray[path-inputIndexBase]=outArray[path-inputIndexBase]+1.0
    return outArray

def get_pathData_node_frequency_array(pathData,nNodes,inputIndexBase=0,verbose=False):
    outArray=get_pathData_node_count_array(pathData,nNodes,inputIndexBase,verbose)/\
        pathData['count']
    return outArray
    

def get_pathData_edge_count_matrix(pathData,nNodes,inputIndexBase=0,verbose=False):
    outMat=np.matrix(np.zeros([nNodes,nNodes]))
    for path in pathData['paths']:
        outMat[path[0:(len(path)-1)]-inputIndexBase,path[1:len(path)]-inputIndexBase]=outMat[
            path[0:(len(path)-1)]-inputIndexBase,path[1:len(path)]-inputIndexBase]+1.0
        outMat[path[1:(len(path))]-inputIndexBase,path[0:(len(path)-1)]-inputIndexBase]=outMat[
            path[1:(len(path))]-inputIndexBase,path[0:(len(path)-1)]-inputIndexBase]+1.0
    return outMat

def get_pathData_edge_frequency_matrix(pathData,nNodes,inputIndexBase=0,verbose=False):
    outMat=get_pathData_edge_count_matrix(pathData,nNodes,inputIndexBase,verbose)/\
        pathData['count']
    return outMat

def serialize_pathData_lengths(pathData):
    return ','.join(map(str,pathData['lengths']))

def serialize_pathData_paths(pathData):
    return ','.join(
        ['_'.join(map(str,pathArray)) for pathArray in pathData['paths']])

def get_1Darray_maxRanks(valArray,unsorted=True,invert=False,verbose=False):
    if unsorted:
        sorting_Array=np.argsort(valArray)
        sortedVals=[valArray[sInd] for sInd in sorting_Array]
        desorting_Array=np.zeros(len(sorting_Array),dtype=int)
        desorting_Array[sorting_Array]=np.arange(len(sorting_Array))
    else:
        sortedVals=valArray
    u,v=np.unique(sortedVals,return_inverse=True)
    if verbose:
        print('u:', end=' ')
        print(u)
        print('v:', end=' ')
        print(v)
    maxRanks=(np.cumsum(np.bincount(v,minlength=u.size))-1)[v]
    if unsorted:
        if invert:
            maxRank=np.max(maxRanks)
            outVals=[maxRank-maxRanks[sInd] for sInd in desorting_Array]
        else:
            outVals=[maxRanks[sInd] for sInd in desorting_Array]
        return outVals
    else:
        if invert:
            maxRank=np.max(maxRanks)
            return [maxRank-maxRanks[sInd] for sInd in np.arange(len(maxRanks))]
        else:
            return maxRanks
        
def get_matrix_element_maxRankings(mat,invert=False,verbose=False):
    unfoldedMat=np.array(copy.deepcopy(mat)).reshape(np.prod(mat.shape))
    unfoldedRankingMat=np.array(get_1Darray_maxRanks(
            unfoldedMat,invert=invert,verbose=verbose))
    rankingMat=unfoldedRankingMat.reshape(mat.shape)
    return rankingMat


