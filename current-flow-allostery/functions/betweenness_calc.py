"""
functions for generating current flow betweenness data from network matrices
"""
import copy
import numpy as np
import scipy as sp
import pandas as pd
import scipy.sparse
import time
import gc

def matLap(mat):
    if mat.shape[0]>mat.shape[1]:
        Lmat=copy.deepcopy(-mat.T).astype(float)
    else:
        Lmat=copy.deepcopy(-mat).astype(float)
    for iRow in np.arange(Lmat.shape[0]):
        Lmat[iRow,iRow]=0
        Lmat[iRow,iRow]=Lmat[iRow,iRow]-\
            np.sum(Lmat[iRow,:])
    if mat.shape[0]>mat.shape[1]:
        Lmat=Lmat.T
    return Lmat

def matAdj(mat):
    if mat.shape[0]>mat.shape[1]:
        Amat=copy.deepcopy(mat.T)
    else:
        Amat=copy.deepcopy(mat)
    Amat[np.arange(Amat.shape[0]),np.arange(Amat.shape[0])]=np.zeros(Amat.shape[0])
    #for iRow in np.arange(Amat.shape[0]):
    #    Amat[iRow,iRow]=0
    if mat.shape[0]>mat.shape[1]:
        Amat=Amat.T
    return Amat
    
def e_btw_from_Linv(Linv,Amat,sources,targets,verbose=False,verboseLevel=0,
                    useLegacyAlgorithm=False,useProgressBar=False):
    if useLegacyAlgorithm:
        return(e_btw_from_Linv_legacy(Linv,Amat,sources,targets,verbose=False,verboseLevel=0,
                    useProgressBar=True))
    else:
        
        if verbose:
            t1=time.time()
        bVecMat=np.zeros((Amat.shape[1],len(sources)*len(targets)))
        smat,tmat=np.meshgrid(sources,targets)
        bVecMat[
            smat.flatten(),np.arange(len(smat.flatten()))
        ]=1
        bVecMat[
            tmat.flatten(),np.arange(len(tmat.flatten()))
        ]=-1
        smat=[]
        tmat=[]
        gc.collect()
        
        potVecMat=np.array(np.matmul(Linv,bVecMat))
        bVecMat=[]
        gc.collect()
        
        nzInds=np.nonzero(Amat)
        btw2=np.array(sp.sparse.coo_matrix(
            (Amat[nzInds]*np.sum(
                np.abs(potVecMat[nzInds[1],:]-potVecMat[nzInds[0],:])/(len(sources)*len(targets)),
                axis=1
            ),
             nzInds),
            shape=Amat.shape
        ).todense())
        if verbose:
            t2=time.time()
            print('total betweenness calculation time:',t2-t1)
        return(btw2)
    
    
def e_btw_from_Linv_legacy(Linv,Amat,sources,targets,verbose=False,verboseLevel=0,
                    useProgressBar=True):
    """
    This algorithm makes use of the full moore-penrose pseudo-inverse
    which is generally quite dense. The algorithm will scale quite poorly
    for very large systems.
    The bottleneck step (after computation of the pseudo inverse) scales as O(m*s*t) where: 
       m is the number of non-zero ntries in the adjacency matrix (network edges)
       s is the number of source nodes
       t is the number of target nodes
    
    This could yield O(n^4) for networks of n nodes in the worst case!
    
    This is only correct for the case where sources and targets are disjoint
    sets. If not, the scaling factor must be adjusted to be equal the number of
    unique combinations of sources and targets.
    """
    eMat=copy.deepcopy(Amat).astype(float)
    #some basic sanity checks
    if((Linv.shape[0]!=Linv.shape[1]) or 
       (Amat.shape[0]!=Amat.shape[1])):
        print("ERROR! Input matrices must by square!")
        eMat[:,:]=0
        return(eMat)
    if((Linv.shape[0]!=Amat.shape[0]) or
       (Linv.shape[1]!=Amat.shape[1])):
        print("ERROR! Input matrices must have the same shape!")
        eMat[:,:]=0
        return(eMat)
    if ((np.min(sources)<0) or
        (np.min(targets)<0) or
        (np.max(sources)>=Linv.shape[0]) or
        (np.max(targets)>=Linv.shape[1])):
        print("ERROR! invalid source or target index detected!")
        eMat[:,:]=0
        return(eMat)
    #get indices of edges... e.g. non-zero entries of Amat
    (Ei,Ej)=np.nonzero(Amat)
    if verbose:
        if verboseLevel > 2:
            print("Linv:", end=' ')
            print(Linv)
            print("Amat:", end=' ')
            print(Amat)
        print("computing betweenness for %g edges"%len(Ei))
        if verboseLevel > 0:
            print("Ei:", end=' ')
            print(Ei)
            print("Ej:", end=' ')
            print(Ej)
    if verbose and useProgressBar:
        Ebtw=np.array(list(map(lambda i,j:
                    Amat[i,j]*\
                    np.sum([np.sum([np.abs(Linv[i,src]+Linv[j,trg]-\
                                   Linv[i,trg]-Linv[j,src]) for trg in targets]) for src in sources]),
                          Ej)))/(len(sources)*len(targets))
    else:
        Ebtw=np.array(list(map(lambda i,j:
                    Amat[i,j]*\
                    np.sum([np.sum([np.abs(Linv[i,src]+Linv[j,trg]-\
                                   Linv[i,trg]-Linv[j,src]) for trg in targets]) for src in sources]),
                          Ei,
                          Ej)))/(len(sources)*len(targets))
    if verbose:
        if verboseLevel > 0:
            print("Ebtw:", end=' ')
            print(Ebtw)
            if verboseLevel > 1:
                print("(Ei,Ej):Ebtw;eMat")
    for iInd in np.arange(len(Ebtw)):
        eMat[Ei[iInd],Ej[iInd]]=Ebtw[iInd]
        if verbose :
            if verboseLevel > 1:
                print("(", end=' ')
                print(Ei[iInd], end=' ')
                print(",", end=' ')
                print(Ej[iInd], end=' ')
                print("):", end=' ')
                print(Ebtw[iInd], end=' ')
                print(";", end=' ')
                print(eMat[Ei[iInd],Ej[iInd]])
    return(eMat)
 
def calcCorrDissipation(corrMat,btwMat):
    return np.sum(np.abs(np.array(btwMat))*\
                  np.abs(np.array(btwMat))/np.abs(np.array(corrMat)))

def getBtwMat(mat,sources,targets,
              verbose=False,verboseLevel=0,
              useProgressBar=False,useLegacyAlgorithm=False):
    """
    Given a (possibly weighted) network in matrix format (mat)
    and a set of source and target nodes (sources and targets)
    return the corresponding network with flow betweenness
    edge weights.
    The pseudo inverse step is the most likely bottleneck, however
    the betweenness calculation iteself scales as O(m*s*t)
    where m=number of network edges, s=number of source nodes, and
    t=number of target nodes. 
    At worst case, this could yield O(n^4) for an n-node matrix!
    Also, the sources and targets must be disjoint sets or the results
    will be incorrect.
    """
    if verbose:
        print("computing matrix Laplacian")
    Lmat=matLap(copy.deepcopy(mat))
    if verbose:
        print("extracting weighted adjacency matrix")
    Amat=matAdj(copy.deepcopy(mat))
    if verbose:
        print("computing moore-penrose inverse of matrix Laplacian")
    Linv=np.linalg.pinv(Lmat)
    if verbose:
        print("generating flow betweenness scores")
    return(e_btw_from_Linv(Linv,Amat,sources,targets,
                           verbose=verbose,verboseLevel=verboseLevel,
                           useLegacyAlgorithm=useLegacyAlgorithm,useProgressBar=useProgressBar))


