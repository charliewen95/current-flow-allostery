#taken directly from the networkx manual
def k_shortest_paths(G, source, target, k, weight=None):
     return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))
    
def converge_subopt_paths_betweenness(inputNetwork,source,target,weight='weight',
                                maxPaths=100,tolerance=1e-6,giveAlphas=False,verbose=False):
    '''Take additional paths between a source / target pair until the betweenness
       centrality of nodes within those paths computed over the attained suboptimal paths
       converges to a given tolerance.
       the convergence critera 'alpha' is computed as:
       alpha = (sum(current path node usage counts) / (number of paths + sum(all path node usage counts))
       the function will return a list of all generated paths. If the "giveAlphas" option is
       turned on, it will also return a list of the alpha after each iteration (path generated)
       this is useful for comparing convergence over different weighting schemes.
       if the "verbose" option is turned on, the alpha after each iteration will be printed
       to standard out as calculation proceeds.'''
    nNodes=len(inputNetwork.nodes())
    pathGenerator=nx.shortest_simple_paths(inputNetwork, source, target, weight=weight)
    nodeCounts=np.zeros(nNodes)
    pathList=[]
    alphas=[]
    newPath=next(pathGenerator)
    pathList.append(copy.deepcopy(newPath))
    pathCounts=np.unique(newPath,return_counts=True)
    tempCounts=np.zeros(nNodes)
    tempCounts[pathCounts[0]]=pathCounts[1]
    nPaths=1
    alpha=np.sum(tempCounts)/np.sum(nodeCounts+tempCounts)/nPaths*1.
    alphas.append(alpha)
    if verbose:
            print("%.3e"%alpha,end=", ")
    while (alpha>tolerance) & (nPaths<maxPaths):
        nodeCounts+=tempCounts
        
        newPath=next(pathGenerator)
        pathList.append(copy.deepcopy(newPath))
        pathCounts=np.unique(newPath,return_counts=True)
        tempCounts=np.zeros(nNodes)
        tempCounts[pathCounts[0]]=pathCounts[1]
        nPaths+=1
        alpha=np.sum(tempCounts)/np.sum(nodeCounts+tempCounts)/nPaths*1.0
        alphas.append(alpha)
        if verbose:
            print("%.3e"%alpha,end=", ")
    if alpha>tolerance:
        print("Maximum number of paths reached before converging betweenness score")
        print("Last relative count delta: %.3e"%alpha)
    if giveAlphas:
        return((pathList,alphas))
    else:
        return(pathList)
 

#A convenience function for calculating the length of an arbitrary path
#in a weighted graph
def calculatePathLength(pathGraph,path,weight='weight'):
    return(np.sum([pathGraph.edges()[(edge[0],edge[1])][weight] \
                   for edge in zip(path[:-1],path[1:])]))


