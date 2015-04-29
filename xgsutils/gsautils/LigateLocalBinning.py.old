#!/usr/bin/env python
import subprocess
import numpy as np

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
def LargebinWrite(prefix,allbins,largebin,id0,K0):
    # build the indicator matrix
    binind = ()
    for i in np.arange(K0):
        for x in largebin[id0][i]:
            binind += (x,)
    
    # count the bin frequency
    binind = np.asarray(binind)
    binpos = sorted(largebin.keys())
    binfreq = dict({pos:{} for pos in binpos})
    for t in np.arange(len(binpos)):
        bc = np.bincount(binind[:,t])
        binfreq[binpos[t]] = 1./(bc+1e-8)

    # iterate over bins
    bc = -1
    for i in np.arange(K0):
        for binvec in largebin[id0][i]:
            bc = bc+1
            binR = LargebinDump(allbins,binpos,binfreq,binvec)
            with open(prefix+'_'+str(binpos[0])+'_'+str(binpos[-1])+'_bin'+str(bc),'w') as f:
                f.write('\n'.join(binR))


from itertools import product
#=======================================
def LigateLocalBin(genome,prefix):
    # find all local binning files
    p = subprocess.Popen(['ls -1 '+'*_'+genome+'_*.long_read_pairs'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    allbins = dict()
    binfiles,err = p.communicate()
    binfiles = binfiles.split()
    for binfile in binfiles:
        id,K,info = ReadLocalBin(binfile)
        allbins[id] = info

    skeys = sorted(allbins.keys())
    # aggregation
    # counting large-bin
    c = 0    # counting large-bin
    # data structure LargeBin
    largebin = dict({skeys[0]:{}})
    # initialize at first position
    largebin[skeys[0]].update({i:([i],) for i in np.arange(allbins[skeys[0]]['K'])})
    # scan through other positions
    for id in skeys[1:]:
        # previous position
        skeys2 = sorted(largebin.keys())
        id0    = skeys2[-1]
        K0     = allbins[id0]['K']
        # current position
        K      = allbins[id]['K']
        # adjacency matrix
        M = np.zeros((K0,K))
        # compute adjacency matrix
        for i,j in product(np.arange(K0),np.arange(K)):
            A = len(allbins[id0]['bins'][i] & allbins[id]['bins'][j])
            B = len(allbins[id0]['bins'][i] | allbins[id]['bins'][j])
            M[i,j] = A/(B+0.)

        # print out if K0<K
        if K0<K and len(largebin)>1:
            # print
            LargebinWrite(prefix,allbins,largebin,id0,K0)
            # build the indicator matrix
            #binind = ()
            #for i in np.arange(K0):
            #    for x in largebin[id0][i]:
            #        binind += (x,)
            # 
            # count the bin frequency
            #binind = np.asarray(binind)
            #binpos = sorted(largebin.keys())
            #binfreq = dict({pos:{} for pos in binpos})     
            #for t in np.arange(len(binpos)):
            #    bc = np.bincount(binind[:,t])
            #    binfreq[binpos[t]] = 1./bc
            #
            # iterate over bins
            #bc = -1
            #for i in np.arange(K0):
            #    for binvec in largebin[id0][i]:
            #        bc = bc+1
            #        binR = LargebinDump(allbins,binpos,binfreq,binvec)
            #        with open(prefix+'_'+str(binpos[0])+'_'+str(binpos[-1])+'_bin'+str(bc),'w') as f:
            #            f.write('\n'.join(binR))
            #
            # reset largebin
            c = c+1
            largebin = dict({id0:{},id:{}})
            largebin[id0].update({i:([i],) for i in np.arange(K0)})
            largebin[id].update({j:() for j in np.arange(K)})
            # column-wise assignmant
            I = np.squeeze(np.asarray(M.argmax(axis=0)))
            J = np.arange(K)
            for i,j in zip(I,J):
                for x in largebin[id0][i]:
                    largebin[id][j] += (x+[j],)

        else:
            # maxiaml assignmnet
            # case 1
            if K0<K:
                # column-wise assignment
                I = np.squeeze(np.asarray(M.argmax(axis=0)))
                J = np.arange(K)
                largebin[id] = dict({j:() for j in np.arange(K)})
                for i,j in zip(I,J):
                    for x in largebin[id0][i]:
                        largebin[id][j] += (x+[j],)

            elif K0>K:
                # row-wise assignment
                I = np.arange(K0)
                J = np.squeeze(np.asarray(M.argmax(axis=1)))
                largebin[id] = dict({j:() for j in np.arange(K)})
                for i,j in zip(I,J):
                    for x in largebin[id0][i]:
                        largebin[id][j] += (x+[j],)

            else:
                # column-wise assignment
                I = np.squeeze(np.asarray(M.argmax(axis=0)))
                J = np.arange(K)
                # collide
                if len(np.unique(I)) < K:
                     # print
                     LargebinWrite(prefix,allbins,largebin,id0,K0)
                     #for i in np.arange(K0):
                     #    for x in largebin[id0][i]:
                     #        print c,': ',','.join(['%s:%s'%(a,b) for a,b in zip(sorted(largebin.keys()),x)])

                     # reset largebin
                     c = c+1
                     largebin = dict({id:{}})
                     largebin[id].update({j:([j],) for j in np.arange(K)})
                # not collide
                else:
                    largebin[id] = dict({j:() for j in np.arange(K)})
                    for i,j in zip(I,J):
                        for x in largebin[id0][i]:
                            largebin[id][j] += (x+[j],)

    if len(largebin) > 0:
        LargebinWrite(prefix,allbins,largebin,id,K)
        #for i in np.arange(K):
        #    for x in largebin[id][i]:
        #        print c,': ',','.join('%s:%s'%(a,b) for a,b in zip(sorted(largebin.keys()),x))



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
