#!/usr/bin/env python

import numpy as np
import networkx as nx
from networkx.algorithms import bipartite
import scipy.linalg as la
from scipy.cluster.vq import vq,kmeans,whiten
from scipy.cluster.hierarchy import fclusterdata
from scipy.spatial.distance import cdist,pdist,squareform
import itertools
from sklearn.cluster import KMeans

#================================
def GraphBuild(dist_mat):
    # get the column number
    with open(dist_mat,'r') as f:
        ncol = len(f.readline().split())

    # read in the edit distance matrix
    B = np.genfromtxt(dist_mat,names=True,usecols=range(1,ncol),dtype=float)
    # the name of local sequences
    Local_Seq = np.array(B.dtype.names)
    # the name of long sequences
    Long_Seq = np.array(np.genfromtxt(dist_mat,skip_header=0,usecols=0,dtype=str))
    # convert ndarray to 2d array
    C = np.array([list(x) for x in B])
    C = C.T
    # remove null local-seq
    K = np.squeeze(np.asarray(C.sum(1)))
    K = K/K.sum()
    I = np.where(K<0.001)
    if len(I) is not 0:
        Local_Seq = np.delete(Local_Seq,I)
        C = np.delete(C,I,0)

    # create bipartite graph
    G = nx.Graph()
    G.add_nodes_from(Local_Seq,bipartite=0)
    G.add_nodes_from(Long_Seq,bipartite=1)
    G.add_edges_from([(Local_Seq[i],Long_Seq[j],{'weight':np.exp(-0.01*(C[i,j]-1)/2.)}) for i,j in zip(*C.nonzero())])
    return (G,Local_Seq,Long_Seq)


#===================================
def GraphSVD(G,row_nodes,column_nodes):
    A = bipartite.biadjacency_matrix(G,row_order=row_nodes,column_order=column_nodes)
    D1 = np.diag(np.squeeze(np.asarray(A.sum(1))))
    D2 = np.diag(np.squeeze(np.asarray(A.sum(0))))
    An = la.sqrtm(la.inv(D1))*A*la.sqrtm(la.inv(D2))
    U,S,V = np.linalg.svd(An)
    V = V.T
    return (U,S,V,D1,D2)


#===================================
def GraphRowNodeClustering0(row_nodes,prefix1,prefix2,K,Z,thresh):
    k = np.arange(K)
    # Z = np.round(la.sqrtm(la.inv(D1))*U[:,k],8)
    # clustering
    C = fclusterdata(Z,thresh)
    # number of cluster
    n = C.max()
    # information of cluster
    C_info = {}
    split = False
    for c in np.arange(n):
        k = np.where(C==(c+1))[0]
        # centroid
        cc = np.mean(Z[k,:],0)
        # list of data
        r1=[]
        r2=[]
        for i in k:
            s = row_nodes[i]
            if s.startswith(prefix1):
                r1.append(i)
            
            if s.startswith(prefix2):
                r2.append(i)
        
        is_split = False;
        if len(r1)>1 and len(r2)!=0:
            is_split = True;

        if len(r1)!=0 and len(r2)>1:
            is_split = True;
 
        split |= is_split
        C_info.update({c+1:{'centroid':cc,'r1':r1,'r2':r2,'split':is_split}})

    return (C_info,split)   


#================================
def GraphRowNodeClustering(row_nodes,prefix1,prefix2,S,U,D1,thresh):
    split = True
    k = np.where(S>=0.97)[0]
    K = len(k)-1
    while split:
        K += 1
        Z = np.round(la.sqrtm(la.inv(D1))*U[:,np.arange(K)],8)
        C_info,split = GraphRowNodeClustering0(row_nodes,prefix1,prefix2,K,Z,thresh)

    K = len(C_info)
    # remove something
    for k in np.arange(K):
        if len(C_info[k+1]['r1'])==0 or len(C_info[k+1]['r2'])==0:
            del C_info[k+1]
    
    return C_info


#================================
def BipartiteGraphClustering(dist_mat,prefix1,prefix2,thresh):
    # build the graph
    G,Local_Seq,Long_Seq = GraphBuild(dist_mat)
    # compute graph laplacian spectrum
    U,S,V,D1,D2 = GraphSVD(G,Local_Seq,Long_Seq)
    # clustering
    C_info = GraphRowNodeClustering(Local_Seq,prefix1,prefix2,S,U,D1,thresh)
    # clustering on long seq
    K = len(C_info)
    # report
    print K


#=================================
def BipartiteGraphClusteringK(dist_mat,prefix1,prefix2,K,thresh):
    # build the graph
    G,Local_Seq,Long_Seq = GraphBuild(dist_mat)
    # compute graph laplacian spectrum
    U,S,V,D1,D2 = GraphSVD(G,Local_Seq,Long_Seq)
    # update K if necessary
    #if ((S>0.7).sum()<K):
    #    K = (S>0.7).sum()
    #
    # data matrix
    Z0 = np.round(la.sqrtm(la.inv(D1))*U[:,np.arange(K)],8)
    Z1 = np.round(la.sqrtm(la.inv(D2))*V[:,np.arange(K)],8)
    # compute the centroid
    estimator = KMeans(init='k-means++', n_clusters=K, n_init=10)
    estimator.fit(Z0)
    B = estimator.cluster_centers_
    # compute the distance between long seq and centroid
    Y = cdist(Z1,B)
    # compute the assignment of long seq
    C = Y.argmin(axis=1)
    # output the pairs of long-seqs
    for k in np.arange(K):
        print "#%d"%k + "\t" + ":".join(['%.2f']*len(S)) % tuple(S)
        I = np.where(C==k)[0]
        for a in Long_Seq[I]:
            print a


from optparse import OptionParser
#================================
def main():
    usage = "Usage: %prog [--t <FLT>] [--k <INT>] --tbl <TABLE> --p1 <PREFIX1> --p2 <PREFIX2> \n\n" \
            "Author: Zeng, Feng; zengfeng@xmu.edu.cn \n" \
            "Last Upate: 2015-03-18 \n\n" \
            "Description: It is to find the connected components in the bipartite graph" 
    parser = OptionParser(usage)
    parser.add_option("--tbl",dest="table",type="string",help="biadjacency-like table",metavar="FILE")
    parser.add_option("--p1",dest="prefix1",type="string",help="prefix name",metavar="STR")
    parser.add_option("--p2",dest="prefix2",type="string",help="prefix name",metavar="STR")
    parser.add_option("--t",dest="thresh",type="float",help="threshold value for k-means",metavar="FLT")
    parser.add_option("--k",dest="num_of_cluster",type="int",help="the number of cluster",metavar="INT")
    (options,args) = parser.parse_args()
    if (options.table is None) or (options.prefix1 is None) or (options.prefix2 is None):
        parser.print_help()
        exit(-1)

    if options.num_of_cluster is None:
        BipartiteGraphClustering(options.table,options.prefix1,options.prefix2,options.thresh)
    else:
        BipartiteGraphClusteringK(options.table,options.prefix1,options.prefix2,options.num_of_cluster,options.thresh)

#================================
if __name__ == "__main__":
    main()
