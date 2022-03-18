"""

"""
from __future__ import absolute_import
from __future__ import print_function

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
import re
import collections
import glob

#Others
import networkx as nx

#self defined functions
#from . import correlation_data_utilities as corr_utils
#from . import database_mod as db_m

def is_number(val):
    try:
        float(val)
    except ValueError:
        return False
    return True

def read_gCorrelation_data(filePath,verbose=False):
    with open(filePath) as inputFile:
        line = inputFile.readline()
        #lineTokens=line.split()
        while line:
            if("[" in line):
                foundStartLine=True
            if foundStartLine:
                lineTokens=line.split()
                if(lineTokens[0].isdigit() & lineTokens[2].isdigit()):
                    nRows=int(lineTokens[0])
                    nCols=int(lineTokens[2])
                break
        line=inputFile.readline()
        nEntries=nRows*nCols
        entries=np.zeros(nEntries)
        count=0
        if(nEntries > 1000):
            cMod=np.ceil(nEntries/100)
            lCount=0
        elif nEntries > 100:
            cMod=np.ceil(nEntries/10)
        else:
            cMod=1
        if verbose:
            print("Reading file "+filePath+":", end=' ')
        while line:
            lineTokens=line.split()
            lineEntries=np.extract(list(map(is_number,np.array(lineTokens))),np.array(lineTokens))
            entries[count:(count+len(lineEntries))]=lineEntries
            count=count+len(lineEntries)
            line=inputFile.readline()
            if verbose:
                if count % cMod == 0:
                    print(str(np.floor(1000 * count / nEntries)/10)+"%", end=' ')
        if verbose:
            print("")
    return collections.OrderedDict({"nRows":nRows,"nCols":nCols,"entries":entries})

def write_carma_matrix(filepath,nRows,nCols,dataEntries,
                       writeHeader=True,
                       header="((protein) and (not hydrogen)) and ((name CA) )",
                       verbose=False):
    if verbose:
        print("Writting carma matrix format file: "+filepath)
    with open(filepath,'w') as outputFile:
        if writeHeader:
            if verbose:
                print('Writting header')
            outputFile.write(header)
            outputFile.write('\n')
        if verbose:
            print('Writting '+str(nRows)+' rows:', end=' ')
            if nRows > 100:
                cMod=nRows/10
            elif nRows > 10000:
                cMod=nRows/100
                count=0
            else:
                cMod=1
        for iRow in np.arange(nRows):
            iStart=iRow*nCols
            iEnd=(iRow+1)*nCols
            outputFile.write(' '.join(map(str,dataEntries[iStart:iEnd])))
            outputFile.write('\n')
            if verbose and (iRow % cMod == 0):
                if cMod>1:
                    if (nRows > 10000) and (count > 9):
                        print("\,    ")
                        count=0
                    print("%4.1f%s "%((np.ceil((1000*iRow)/nRows)/10.0),"%"), end=' ')
                    if nRows > 10000:
                        count = count+1
                else:
                    print(".", end=' ')
        if verbose:
            print("\ndone")
            
            

def read_carma_matrix(filepath,has_header=True,returnHeader=False,verbose=False):
    #assumes all carma matrices are square matrices!
    if verbose:
        print("Reading carma matrix format file: "+filepath)
    with open(filepath,'r') as inputFile:
        if(has_header):
            if verbose:
                print('reading header')
            headerLine=inputFile.readline()
        line=inputFile.readline()
        count=0
        cMod=1
        while line:
            lineTokens=line.split()
            if(count==0):
                nRows=len(lineTokens)
                nEntries=nRows*nRows
                entries=np.zeros(nEntries)
                if(verbose):
                    print('Reading '+str(nRows)+' rows: ', end=' ')
                    if nRows > 1000:
                        cMod = np.ceil(nRows / 100)
                    elif nRows > 100:
                        cMod = np.ceil(nRows / 10)
                    else:
                        cMod = 1
            valArray=np.array(lineTokens)
            if (len(valArray)==nRows):
                iStart=count*nRows
                iEnd=iStart+len(valArray)
                sys.stdout.flush()
                entries[iStart:iEnd]=valArray
                if verbose & (count % cMod == 0):
                    if nRows > 10:
                        print(str(np.floor(1000*count/nRows)/10)+"% ", end=' ')
                    else:
                        print(".", end=' ')
                count=count+1
            line=inputFile.readline()
        if verbose:
            print("")
    if count > 0:
        if has_header & returnHeader:
            return collections.OrderedDict({"headerLine":headerLine,"nRows":nRows,"nCols":nRows,
                                            "entries":entries})
        else:
            return collections.OrderedDict({"nRows":nRows,"nCols":nRows,"entries":entries})
    else:
        print("Error! Data file appears empty!")

def convert_gCorrelationData_to_carmaMatrix(gCorrFilePath,outputFilePath,
                                            writeHeader=True,
                                            header="((protein) and (not hydrogen)) and ((name CA) )",
                                            verbose=False):
    #if verbose:
    #    print "Reading g_correlation data file"
    #    sys.stdout.flush()
    gCorrData = read_gCorrelation_data(filePath=gCorrFilePath,verbose=verbose)
    #if verbose:
    #    print "Writting g_correlation data to carma matrix format file"
    #    sys.stdout.flush()
    write_carma_matrix(filepath=outputFilePath,
                       nRows=gCorrData['nRows'],nCols=gCorrData['nCols'],
                       dataEntries=gCorrData['entries'],
                       writeHeader=writeHeader,header=header,verbose=verbose)
    if verbose:
        print("Conversion complete")

   
#A convenience function for calculating the length of an arbitrary path
#in a weighted graph
def calculatePathLength(pathGraph,path,weight='weight'):
    return(np.sum([pathGraph.edges()[(edge[0],edge[1])][weight] \
                   for edge in zip(path[:-1],path[1:])]))

#Utilities for computing distance topology in pytraj
def checkCollision(Rmin1,Rmax1,Rmin2,Rmax2,collisionRadius=0,axis=None):
    return(np.product((Rmin1-collisionRadius)<=(Rmax2+collisionRadius),axis=axis) * \
           np.product((Rmax1+collisionRadius)>=(Rmin2-collisionRadius),axis=axis))

def collisionCount(Rmin1,Rmax1,Rmin2,Rmax2,
                   collisionRadius=0,crdAxis=1):
    return(np.sum(np.product((Rmin1-collisionRadius)<=(Rmax2+collisionRadius),axis=crdAxis) * \
                  np.product((Rmax1+collisionRadius)>=(Rmin2-collisionRadius),axis=crdAxis)))

#This just returns whether or not a collision has ever happened, only slightly faster
#than the series apparently, though likely smaller memory requirements
def compute_BoxCollision_matrix(traj,collisionRadius=0.,resInds=None,showProgress=False):
    if resInds is None:
        nRes=traj.top.n_residues()
        resnums=np.arange(nRes)
    else:
        resnums=resInds
        nRes=len(resInds)
    
    resInds=[traj.topology.atom_indices(':%g'%iRes) for iRes in resnums+1]
    
    if showProgress:
        resIter='Computing Residue Minimum Bounds'
    else:
        resIter=resInds
    resMinBounds=np.array([np.min(traj.xyz[:,resInd,:],axis=1) \
                           for resInd in resIter])
    
    if showProgress:
        resIter='Computing Residue Maximum Bounds'
    else:
        resIter=resInds
    resMaxBounds=np.array([np.max(traj.xyz[:,resInd,:],axis=1) \
                           for resInd in resIter])
    
    resPairs=np.array(list(combinations(np.arange(nRes),2)))
    if showProgress:
        pairIter='Computing Collisions'
    else:
        pairIter=resPairs
    collisionCheckArray=[
        checkCollision(resMinBounds[resPair[0]],resMaxBounds[resPair[0]],
                       resMinBounds[resPair[1]],resMaxBounds[resPair[1]],
                       collisionRadius=collisionRadius) \
        for resPair in pairIter]

    collisionMat=sp.sparse.coo_matrix(
        (collisionCheckArray,
         (resPairs[:,0],resPairs[:,1])),shape=(nRes,nRes))
    collisionMat=collisionMat+collisionMat.T
    return(collisionMat)

#Counts the number of frames where each residue pair has collided
#returns the result as a sparse matrix (scipy coo format)
def compute_BoxCollision_CountMatrix(traj,collisionRadius=0.,
                                     resinds=None,
                                     minBounds=None,maxBounds=None,
                                     showProgress=False,
                                     frameAxis=0,indAxis=1,crdAxis=2,
                                     returnBoundVecs=False):
    if resinds is None:
        nRes=traj.top.n_residues
        resnums=np.arange(nRes)
        resInds=[traj.topology.atom_indices(':%g'%iRes) for iRes in resnums+1]
    else:
        #resnums=resInds
        nRes=len(resinds)
        resInds=resinds
    
    #print(len(resInds))
    if (minBounds is None): #| (len(minBounds) != len(resInds)):
        if showProgress:
            resIter='Computing Residue Minimum Bounds'
        else:
            resIter=resInds
        resMinBounds=np.array([np.min(traj.xyz[:,resIndSet,:],axis=indAxis) \
                               for resIndSet in resIter])
    else:
        resMinBounds=minBounds
    
    if (maxBounds is None): #| (len(maxBounds) != len(resInds)):
        if showProgress:
            resIter='Computing Residue Maximum Bounds'
        else:
            resIter=resInds
        resMaxBounds=np.array([np.max(traj.xyz[:,resIndSet,:],axis=indAxis) \
                               for resIndSet in resIter])
    else:
        resMaxBounds=maxBounds
    
    resPairs=np.array(list(combinations(np.arange(nRes),2)))
    if showProgress:
        pairIter='Computing Collisions'
    else:
        pairIter=resPairs
    #print(resMinBounds.shape)
    collisionCheckArray=[
        collisionCount(resMinBounds[resPair[0]],resMaxBounds[resPair[0]],
                       resMinBounds[resPair[1]],resMaxBounds[resPair[1]],
                       collisionRadius=collisionRadius,crdAxis=crdAxis-1) \
        for resPair in pairIter]

    #return(collisionCheckArray)
    collisionMat=sp.sparse.coo_matrix(
        (np.concatenate([collisionCheckArray,collisionCheckArray]),
         (np.concatenate([resPairs[:,0],resPairs[:,1]]),
          np.concatenate([resPairs[:,1],resPairs[:,0]]))),shape=(nRes,nRes))
    if returnBoundVecs:
        return(collisionMat,resMinBounds,resMaxBounds)
    else:
        return(collisionMat)
    
#Utilities to compute pearson and linear mutual information correlation
#matrices... these seem to be slow compared to carma and g_corr but may
#be useful when such tools cannot be easily compiled (e.g. cloud computing applications)
def calc_Ci(X,crdAxis=1):
    return(
        np.mean(
            np.apply_along_axis(
                lambda x: np.array(np.matrix([x]).T*np.matrix([x])),
                arr=X,
                axis=crdAxis),axis=0))

def calc_Cij(Xi,Xj,crdAxis=1):
    return(
        np.mean(
            np.apply_along_axis(
                lambda x: np.array(np.matrix([x]).T*np.matrix([x])),
                arr=np.concatenate([Xi,Xj],axis=crdAxis),
                axis=crdAxis),axis=0))

def calc_Linear_Mutual_Information(Xi,Xj,
                                   Ci=None,Cj=None,
                                   #Cii=None,Cjj=None,
                                   crdAxis=1,
                                   verbose=False):
    #Ci,Cii,Cj, and Cjj can be input if they have been calculated in advance
    #This can save significant when calcuting linear MI over all pairs in a 
    #large number of coordinate sets since Ci,Cii,Cj, and Cjj can be computed
    #in a single loop over all coordinate sets instead of needing to recalculated
    #for each ij coordinate set pair
    if Ci is None:
        ci=calc_Ci(Xi,crdAxis)
    else:
        ci=Ci
    #if Cii is None:
    #    cii=calc_Cij(Xi,Xi,crdAxis)
    #else:
    #    cii=Cii
    if Cj is None:
        cj=calc_Ci(Xj,crdAxis)
    else:
        cj=Cj
    #if Cjj is None:
    #    cjj=calc_Cij(Xj,Xj,crdAxis)
    #else:
    #    cjj=Cjj
        
    cij=calc_Cij(Xi,Xj,crdAxis)
    #cji=calc_Cij(Xj,Xi,crdAxis)
    CijMat=cij
    #CijMat=np.matrix(
    #    np.concatenate(
    #        [np.concatenate([cii,cij],axis=1),
    #         np.concatenate([cij,cjj],axis=1)],
    #        axis=0))
    if verbose:
        for entryName,entry in [
            ['Ci',ci], #['Cii',cii],
            ['Cj',cj], #['Cjj',cjj],
            #['Cij',cij],['Cji',cji],
            ['CijMat',CijMat],
            ['det(Ci)',np.linalg.det(ci)],
            ['det(Cj)',np.linalg.det(cj)],
            ['det(CijMat)',np.linalg.det(CijMat)]
        ]:
            print(entryName)
            print(entry)
    return(.5*(np.log(np.linalg.det(ci))+np.log(np.linalg.det(cj))-np.log(np.abs(np.linalg.det(CijMat)))))

def calc_pear_corr(Xi,Xj,Rii=None,Rjj=None,crdAxis=1,verbose=False):
    if Rii is None:
        rii=np.mean(np.apply_along_axis(lambda x: np.sum(x),arr=Xi**2,axis=crdAxis))
    else:
        rii=Rii
    if Rjj is None:
        rjj=np.mean(np.apply_along_axis(lambda x: np.sum(x),arr=Xj**2,axis=crdAxis))
    else:
        rjj=Rjj
    rij=np.mean(np.apply_along_axis(lambda x: np.sum(x),arr=Xi*Xj,axis=crdAxis))
    if verbose:
        print("rii:",rii)
        print("rjj:",rjj)
        print("rij:",rij)
    return(rij/(np.sqrt(rii)*np.sqrt(rjj)))
    

def calc_Normalized_LinearMI(Xi,Xj,
                             Ci=None,Cj=None,
                             #Cii=None,Cjj=None,
                             Rii=None,Rjj=None,
                             crdAxis=1,verbose=False):
    Imi=calc_Linear_Mutual_Information(Xi,Xj,
                                       Ci,Cj,
                                       #Cii,Cjj,
                                       crdAxis,
                                       verbose)
    Rmi=(1-np.exp(-2*Imi/Xi.shape[crdAxis]))**(1/2)
    
    rij=calc_pear_corr(Xi,Xj,Rii,Rjj,crdAxis)
    Igauss=-Xi.shape[crdAxis]/2.*np.log(1-rij**2)
    Rgauss=(1-np.exp(-2*Igauss/Xi.shape[crdAxis]))**(1/2)
    if verbose:
        print('rij',rij)
        print('Igauss',Igauss)
        print('Rgauss',Rgauss)
        print('Imi',Imi)
        print('Rmi',Rmi)
    return(Rmi)


