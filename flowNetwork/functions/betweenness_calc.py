"""
Functions
---------
.. autofunction:: getBtwMat
"""

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
 
