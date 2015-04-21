#!/bin/python

import argparse
import sys
from scipy import *
from scipy.sparse import *
import numpy as np

def quintal_to_decimal(seq):
    alpha={'A':0,'C':1,'G':2,'T':3,'a':0,'c':1,'g':2,'t':3}
    b=1
    d=0
    for i,k in enumerate(seq):
        d += alpha[k]*b
        if i>0:
            b *= 4

    return d


def kmer_dict(file,ksize,normalize):
    kd = {}
    lines = [line.strip() for line in open(file)]
    z = 0
    for line in lines:
        seq,count = line.split()
        count=float(count)
        if seq in kd:
            kd[seq] += count
        else:
            kd[seq] = count
        z += count

    if normalize == True:
        for seq in kd:
            kd[seq] /= z

    return kd


def similarity(kd1,kd2):
    kmer1 = set(kd1.keys())
    kmer2 = set(kd2.keys())
    kmers = kmer1 | kmer2
    d1=0
    d2=0
    d=0
    for kmer in kmers:
        occ=0
        if kmer in kmer1:
            d1 += kd1[kmer]**2
            occ+=1

        if kmer in kmer2:
            d2 += kd2[kmer]**2
            occ+=1

        if occ==2:
            d += kd1[kmer]*kd2[kmer]

    sim = d/np.sqrt(d1*d2)
    return sim



def seq_cosine_similarity(file1,file2,ksize,normalize):
    kd1 = kmer_dict(file1,ksize,normalize)
    kd2 = kmer_dict(file2,ksize,normalize)
    print similarity(kd1,kd2)

    
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("kmer1",type=str,help="the kmer spectrum of sequence one")
    parser.add_argument("kmer2",type=str,help="the kmer spectrum of sequence two")
    parser.add_argument("-k",dest="ksize",type=int,nargs=1,help="the size of the kmer")
    parser.add_argument("--no-normalize",dest="normalize",action="store_false",help="do not normalize the elements in the vector")
    parser.set_defaults(normalize=True)
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    seq_cosine_similarity(args.kmer1,args.kmer2,args.ksize,args.normalize)
