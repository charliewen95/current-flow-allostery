#Utilities for computing 'flow betweenness' scores for matrix
#representations of networks
def corrDataDictToMat(dataDict):
    return np.matrix(
        np.array(
            dataDict['entries']).reshape(dataDict['nRows'],dataDict['nCols']))

def corrMatToDataDict(mat):
    nRows=mat.shape[0]
    nCols=mat.shape[1]
    entries=np.array(mat).reshape(nRows*nCols)
    return(collections.OrderedDict({
                'nRows':nRows,
                'nCols':nCols,
                'entries':entries}))

def edgeDataDictToMatrix(dataDict,nRows=-1,nCols=-1):
    if nRows < 0:
        nRows=np.max(dataDict['Ei'])+1
    if nCols < 0:
        nCols=np.max(dataDict['Ej'])+1
    outMat=np.matrix(np.zeros([nRows,nCols])*0.0)
    outMat[dataDict['Ei'],dataDict['Ej']]=dataDict['Ew']
    return outMat

def write_dataDict_to_carma_matrix(filepath,dataDict,writeHeader=True,
                                       useDictHeader=True,
                                       header="((protein) and (not hydrogen)) and ((name CA) )",
                                       verbose=False):
    if not ( ('nRows' in dataDict) and \
             ('nCols' in dataDict) and \
             ('entries' in dataDict) ):
        print("ERROR! data dictionary is missing needed key value pairs. !!Aborting!!")
        return
    if (writeHeader and useDictHeader):
        if not ('headerLine' in dataDict) :
            print("WARNING! headerLine was missing from data dictionary.")
            print(" -Defaulting to: "+header)
            headerLine=header
        else:
            headerLine=dataDict['headerLine']
    write_carma_matrix(filepath=filepath,
                       nRows=dataDict['nRows'],nCols=dataDict['nCols'],
                       dataEntries=dataDict['entries'],writeHeader=writeHeader,
                       header=headerLine,verbose=verbose)
        
            

def write_mat_to_carma_matrix(filepath,mat,writeHeader=True,
                              header="((protein) and (not hydrogen)) and ((name CA) )",
                              verbose=False):
    dataDict=corrMatToDataDict(mat)
    if writeHeader:
        dataDict['headerLine']=header
    if writeHeader:
        write_dataDict_to_carma_matrix(filepath,dataDict,
                                       writeHeader=writeHeader,useDictHeader=writeHeader,
                                       header=header,verbose=verbose)
    else:
        write_dataDict_to_carma_matrix(filepath,dataDict,
                                       writeHeader=writeHeader,useDictHeader=writeHeader,
                                       verbose=verbose)


