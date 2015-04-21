#!/bin/python
import sys
import argparse
import numpy as np
import scipy.cluster.hierarchy as hc
from collections import defaultdict

def hcluster(file,t):
    X = []
    Y = []
    lines = [ line.strip() for line in open(file) ]
    for line in lines:
        x1,x2,s = line.split()
        dxy = float(s)
        if x1 not in X:
            X.append(x1)
      
        if x2 not in X:
            X.append(x2)
        
        Y.append(dxy)

    Y = np.array(Y)
    Z = hc.linkage(Y,method='average')
    T = hc.fcluster(Z,t,criterion="distance")
    R = defaultdict(list)    
    for i,t in enumerate(T):
        R[t].append(X[i])
    
    for key,cluster in sorted(R.items(), key=lambda(k,v):len(v), reverse=True):
        print key,' '.join(cluster)


def hcluster2(file,t,hit):
    th=0
    fh={}
    lines=[line.strip() for line in open(hit)]
    for line in lines:
        f,h = line.split()
        fh[f] = int(h)
        th += int(h)

    X = []
    Y = []
    lines = [ line.strip() for line in open(file) ]
    for line in lines:
        x1,x2,s = line.split()
        dxy = float(s)
        if x1 not in X:
            X.append(x1)

        if x2 not in X:
            X.append(x2)

        Y.append(dxy)

    Y = np.array(Y)
    Z = hc.linkage(Y,method='average')
    T = hc.fcluster(Z,t,criterion="distance")
    R = defaultdict(list)
    C = {}
    for i,t in enumerate(T):
        R[t].append(X[i])
        if t in C:
            C[t] += fh[X[i]]
        else:
            C[t] = fh[X[i]]

    qh=0.0
    for key,h in sorted(C.items(), key=lambda(k,v):v, reverse=True):
        qh+=float(h)
        print key,'\t',h,'\t',qh/th,'\t'.join(R[key])


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("distfile",type=str,help="the distance file")
    parser.add_argument("-t",type=float,nargs="?",const=0.01,default=0.01,help="the distance threshold")
    parser.add_argument("--read-hit",type=str,nargs="?",dest="read_hit",help="sorting the clusters by the read hits")

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args()
    df = args.distfile
    t = args.t
    if args.read_hit is None:
        hcluster(df,t)    
    else:
        hcluster2(df,t,args.read_hit)
