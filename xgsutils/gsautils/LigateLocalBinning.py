#!/usr/bin/env python
import subprocess
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

#=======================================
def ReadLocalBin(binfile):
    a = binfile.replace('_',' ').replace('.',' ').split()
    local_name = binfile
    local_p0   = int(a[2])
    local_p1   = int(a[3])
    local_bin  = dict()
    bin_id     = ""
    bin_R      = list()
    K          = 0
    with open(binfile) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                K = K+1
                if len(bin_id) is not 0:
                    local_bin.update({int(bin_id):set(bin_R)})
    
                bin_id = line[1:].split()[0]
                bin_R  = list()
            else:
                bin_R.append(line)

        if len(bin_id) is not 0:
            local_bin.update({int(bin_id):set(bin_R)})

    return (local_p0,K,dict({'name':local_name,'K':K,'p0':local_p0,'p1':local_p1,'bins':local_bin}))


#=======================================
def LargebinDump(allbins,binpos,binfreq,binvec):
    R = set()
    for i in np.arange(len(binpos)):
        N = len(allbins[binpos[i]]['bins'][binvec[i]])
        S = np.random.choice(2,N,p=[1-binfreq[binpos[i]][binvec[i]],binfreq[binpos[i]][binvec[i]]])
        R = R | set(np.asarray(list(allbins[binpos[i]]['bins'][binvec[i]]))[np.where(S==1)[0]])

    return list(R)


#=======================================
#def LargebinWrite(prefix,allbins,largebin,id0,K0):
#    # build the indicator matrix
#    binind = ()
#    for i in np.arange(K0):
#        for x in largebin[id0][i]:
#            binind += (x,)
#    
#    # count the bin frequency
#    binind = np.asarray(binind)
#    binpos = sorted(largebin.keys())
#    binfreq = dict({pos:{} for pos in binpos})
#    for t in np.arange(len(binpos)):
#        bc = np.bincount(binind[:,t])
#        binfreq[binpos[t]] = 1./(bc+1e-8)
#
#    # iterate over bins
#    bc = -1
#    for i in np.arange(K0):
#        for binvec in largebin[id0][i]:
#            bc = bc+1
#            binR = LargebinDump(allbins,binpos,binfreq,binvec)
#            with open(prefix+'_'+str(binpos[0])+'_'+str(binpos[-1])+'_bin'+str(bc),'w') as f:
#                f.write('\n'.join(binR))


#======================================
def LargebinWrite(prefix,allbins,G,L,c):
    R = set()
    for l in L:
        b = G.node[l]['k']
        i = G.node[l]['i']
        R = R | set(np.asarray(list(allbins[b]['bins'][i])))
    
    with open(prefix+"_largebin_"+str(c),'w') as f:
        f.write('\n'.join(R))
    

#======================================
def IsRepairableCollision(A,n):
    C = np.bincount(A)
    return np.size(np.where((C>1) & (C<n)))>0


#======================================
def IsRepairableBreak(A,n):
    C = np.bincount(A)
    return (np.where(C==0)[0],np.size(np.where(C==0))==1)
    

from itertools import product
#=======================================
def LigateLocalBin(genome,prefix):
    ##### create a graph
    G = nx.Graph()
    ##### find all local binning files
    p = subprocess.Popen(['ls -1 '+'*_'+genome+'_*.long_read_pairs'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    allbins = dict()
    binfiles,err = p.communicate()
    binfiles = binfiles.split()
    for binfile in binfiles:
        id,K,info = ReadLocalBin(binfile)
        allbins[id] = info

    skeys = sorted(allbins.keys())
    ##### aggregation
    # counting large-bin
    c = 0    # counting large-bin
    # data structure LargeBin
    largebin = dict({skeys[0]:{}})
    ### initialize at first position
    largebin[skeys[0]].update({i:([i],) for i in np.arange(allbins[skeys[0]]['K'])})
    # add graph node
    for i in np.arange(allbins[skeys[0]]['K']):
        G.add_node(str(skeys[0])+"_"+str(i),{'k':skeys[0],'i':i})

    id0 = skeys[0]
    ### scan through other positions
    for ii in np.arange(1,len(skeys)):
    #for id in skeys[1:]:
        # previous position
        id0    = skeys[ii-1]
        K0     = allbins[id0]['K']
        # current position
        id     = skeys[ii]
        K      = allbins[id]['K']
        # add node at current position
        for i in np.arange(K):
            G.add_node(str(id)+"_"+str(i),{'k':id,'i':i})

        # adjacency matrix
        M = np.zeros((K0,K))
        # compute adjacency matrix
        for i,j in product(np.arange(K0),np.arange(K)):
            A = len(allbins[id0]['bins'][i] & allbins[id]['bins'][j])
            B = len(allbins[id0]['bins'][i] | allbins[id]['bins'][j])
            if A/(B+0.)>=0.01:
                M[i,j] = A/(B+0.)

        # column-wise assignmant
        I = np.squeeze(np.asarray(M.argmax(axis=0)))
        J = np.arange(K)
        II = np.array(I)
        for i,j in zip(I,J):
            if M[i,j] < 0.01:
                II[j] = K0

        # case 1: break
        BI,IB = IsRepairableBreak(II,K0)
        if IB == True:
            I[np.where(II==K0)[0]] = BI
            M[BI,np.where(II==K0)[0]] = 1.0

        # case 2: collide
        if IsRepairableCollision(II,K0):
            if np.size(np.where(II==K0))>0:
                I[np.where(II==K0)] = -1

            for xi in np.arange(K0):
                if np.size(np.where(II==xi))>1:
                    MM = M[xi,np.where(II==xi)[0]].max()
                    for xj in np.where(II==xi)[0]:
                        if M[xi,xj] < MM:
                            I[xj] = -1

        for i,j in zip(I,J):
            if i != -1:
                if M[i,j] >= 0.01:
                    G.add_edge(str(id0)+"_"+str(i),str(id)+"_"+str(j))

    ##### dump connected components
    c = 0
    for L in nx.connected_components(G):
        if len(L)>1:
            LargebinWrite(prefix,allbins,G,L,c)
            c = c+1
   
    #if len(largebin) > 0:
    #    LargebinWrite(prefix,allbins,largebin,id,K)
    #    #for i in np.arange(K):
    #    #    for x in largebin[id][i]:
    #    #        print c,': ',','.join('%s:%s'%(a,b) for a,b in zip(sorted(largebin.keys()),x))



from optparse import OptionParser
#=======================================
def main():
    usage = "Usage: %prog <GENOME_NAME> <PREFIX>\n\n" \
            "Author: Zeng, Feng; zengfeng@xmu.edu.cn\n" \
            "Last Update: 2015-03-31\n\n" \
            "Description: It is to ligate the logcal binning of the given genome"
    parser = OptionParser(usage)
    options,args = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        exit(-1)

    LigateLocalBin(args[0],args[1])

#=======================================
if __name__ == "__main__":
    main()
