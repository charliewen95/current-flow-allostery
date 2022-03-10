#functions for generating network topology and edge weighting data with pytraj
    
def gen_perResCOM_traj(inputTrajPath,inputTopPath,residList,
                              resSelectMask='',resAtomMask='@CA',
                              COM_mask='',computeCommand='vector center :',
                              threads=1,verbose=False):
    #if noPytrajFlag:
    #    print "ERROR! gen_perResCOM_traj cannot be run if pytraj fails to load. Aborting"
    #    return
    #Default is to return per-residue center of mass mapped to alpha carbons of each residue
    if verbose:
        print("loading input trajectory")
    iterTraj=pt.iterload(inputTrajPath,top=inputTopPath)
    if verbose:
        print(iterTraj)
    commandList=["vector center :"+str(rid)+COM_mask for rid in residList]
    if verbose:
        if len(commandList) >= 10:
            print("first 10 mapped commands: ", end=' ')
            print(commandList[0:10])
        else:
            print("mapped commands: ", end=' ')
            print(commandList)
        print("running mapped commands")
    perResData=np.array(list(pt.compute(commandList,iterTraj,n_cores=threads).values()))
    if resSelectMask=='':
        residString=[str(rid) for rid in residList]
        resSelectMask=':'+','.join(residString)+resAtomMask
    if verbose:
        print("extracting residue selection subtrajectory")
    tempTraj=iterTraj[0:iterTraj.n_frames,resSelectMask]
    if verbose:
        print("subtrajectory info:")
        print(tempTraj)
    if not (tempTraj.shape[0]==perResData.shape[1] and \
            tempTraj.shape[1]==perResData.shape[0]):
        print("ERROR! mismatch between trajectory subselection coordinates and perResidue command data")
        print(" trajectory.xyz.shape=", end=' ')
        print(trajectory.xyz.shape, end=' ')
        print("; perResData.shape=", end=' ')
        print(perResData.shape[[1,0,2]])
        return
    if verbose:
        print("updating trajectory coordinates")
    for iDim in np.arange(3):
        tempTraj.xyz[:,:,iDim]=perResData[:,:,iDim].T
    print("done")
    return tempTraj

def gen_smoothed_contact_map(traj,resids,
                             timeAggFun=lambda timeSeries: 1.0*(np.mean(timeSeries)>.75),
                             distSmoothFun=lambda distArray:(6.0-np.clip(distArray,3.0,6.0))/(6.0-3.0),
                             verbose=False,verboseLevel=0):
    #if noPytrajFlag:
    #    print "ERROR! gen_perResCOM_traj cannot be run if pytraj fails to load. Aborting"
    #    return
    #traj should by a pytraj.trajectory or trajectory iterator
    #resids should be a list of residue ids (using AMBER numbering)
    #distSmoothFun should apply a smoothing kernel over a
    #len(resids) by traj.n_frames 2D np.ndarray
    #It defaults to a linear dropoff from an upper cutoff of 3.0 (yields 1.0 or full contact)
    #to a lower cutoff of 6.0 (yields 0.0 or no contact)
    #timeAggFun should apply an appropriate aggregation function
    #it defaults to computing the mean value and using a .75 high pass filter
    #to recover the unsmoothed contact map use distSmoothFun=lambda distArray: 1.0*(distArray<3.0)
    #verbose will turn on or off a simple text based progress monitor
    nResidues=len(resids)
    tempMat=np.zeros(nResidues)
    residArray=np.array(resids)
    if verbose:
        print("Computing contact matrix:", end=' ')
        if nResidues > 100:
            cMod=nResidues/100
        elif nResidues > 10:
            cMod=nResidues/10
        else:
            cMod=1
        lCount=0
    (mgrid_i,mgrid_j)=np.mgrid[0:nResidues,0:nResidues]
    tempMat=np.zeros([nResidues,nResidues])
    for iRow in np.arange(nResidues):
        if verbose:
            if((iRow % cMod == 0) or (iRow==(nResidues-1))):
                print("%4.1f%s "%(np.ceil(1000*iRow/nResidues)/10,"%"), end=' ')
                if lCount==10:
                    print("\n                         ", end=' ')
                    lCount=0
                else:
                    lCount=lCount+1
        rgrid_i=residArray[mgrid_i[iRow,:].flatten()]
        rgrid_j=residArray[mgrid_j[iRow,:].flatten()]
        commandList=list(map(lambda rid,rjd: 'nativecontacts :'+\
                        str(rid)+' :'+str(rjd)+' mindist',rgrid_i,rgrid_j))
        tempDists=np.array(list(pt.compute(commandList,traj).values())[2:(3*len(rgrid_i)):3])
        if verbose and (verboseLevel>0):
            print(commandList[0:3])
            print(type(tempDists), end=' ')
            print(" ", end=' ')
            print(tempDists.shape, end=' ')
            print(tempDists)
        matVals=[timeAggFun(timeSeries) for timeSeries in distSmoothFun(tempDists)]
        tempMat[mgrid_i[iRow,:].flatten(),mgrid_j[iRow,:].flatten()]=matVals
        tempMat[iRow,iRow]=0
    print("")
    return tempMat

def compute_pairwise_minDist_data(traj,
                                  resindexpairs=None,
                                  chunkSize=1000,
                                  outFilePath=None,
                                  returnData=True,
                                  showProgress=False):
    '''
        traj: a pytraj trajectory
        resindexpairs: pairs of residues to compute distance series for. These should be base 0 indices.
            The default will compute all possible residue pair combinations. This option can allow you
            to filter which pairs to compute. E.g. if you have an nResidue by nResidue matrix (filterMatrix)
            wich has nonzero values for only the pairs to be computed you could use:
            resindexpairs=zip(np.nonzero(filterMatrix)[0],np.nonzero(filterMatrix)[1])
        showProgress: if set to True, display a tqdm_notebook style progress bar
        chunkSize: Due to technical considerations, it is apparently faster to do several hundred to
            several thousand pair distance calculations at a time. This controls the number of
            pair distances calculated in one pass over the entire trajectory. For the IGPS system
            1000 seemed to be roughly optimal, but the optimal value will likely vary depending on system size
            and trajectory length.
        outFilePath: if set, the data will be written chunk by chunk to the specified filepath as it is
            generated.
        returnData: return an array containing the computed data. If turned off, nothing will be returned.
            This can be useful when handling a very large number of pairs (i.e. if the computed data
            will not fit in memory). Setting this to false and providing an outFilePath will cause the 
            data to be written directly to disk instead of stored in an in-memory array.
            Note that you will still need to be able to fit at least pair set chunk worth of data in memory.
        
        returns: an M-Residue_Pair by N-Frame array where rows are the residue pair and columns are 
            trajectory frames. Each entry is the minimum residue-residue interatomic distance for the
            given residue pair at the given frame.
    '''
    if resindexpairs is None:
        resIndices=np.arange(traj.Top.n_residues)
        resIndexPairs=list(combinations(resIndices,2))
    else:
        resIndices=resindexpairs
    
    distData=[]
    
    chunkStart=0
    nPairs=len(resIndices)
    
    if showProgress:
        pbar='computing pair distances'
    count=0
    while chunkStart<nPairs:
        chunkEnd=chunkStart+chunkSize
        pairIter=resIndices[chunkStart:chunkEnd]
        commandList=['nativecontacts mindist :{:g} :{:g}'.format(resPair[0]+1,resPair[1]+1) \
                     for resPair in pairIter]
        tempData=list(
            pt.compute(commandList,traj).values())[2::3]
        if returnData:
            distData.append(copy.deepcopy(tempData))
        chunkStart=chunkEnd
        if not (outFilePath is None):
            pbar.set_description('writting data to disk')
            if count==0:
                with open(outFilePath,'w') as outFile:
                    np.savetxt(outFile,X=tempData)
                outFile=open(outFilePath,'a')
            else:
                np.savetxt(outFile,X=tempData)
            pbar.set_description('computing pair distances')
        count+=1
        gc.collect()
        pbar.update(chunkSize)
    pbar.close()
    outFile.close()
    
    if returnData:
        distData=np.concatenate(distData,axis=0)
        return(distData)

def drawProtNetEdge(protStruc,resID1,resID2,ngViewOb,
                    frame=0,edgeColor=[.5,.5,.5],radius=1,
                    *shapeArgs,**shapeKwargs):
    crd1=pt.center_of_mass(protStruc,':%g@CA'%resID1)[frame]
    crd2=pt.center_of_mass(protStruc,':%g@CA'%resID2)[frame]
    
    
    resname1=protStruc.topology.residue(resID1-1).name
    resid1=protStruc.topology.residue(resID1-1).original_resid
    
    resname2=protStruc.topology.residue(resID2-1).name
    resid2=protStruc.topology.residue(resID2-1).original_resid
    edgeLabel='%s.%g-%s.%g (%g-%g)'%(
        resname1,resid1,resname2,resid2,
        resID1-1,resID2-1)
    
    return ngViewOb.shape.add_cylinder(
        list(crd1),list(crd2),edgeColor,radius,
        edgeLabel,
        *shapeArgs,**shapeKwargs)


