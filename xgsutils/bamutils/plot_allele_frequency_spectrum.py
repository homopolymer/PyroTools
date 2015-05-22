#!/usr/bin/env python

from optparse import OptionParser
import sys

import pysam
import pysamstats
import matplotlib.pyplot as plt
import pandas as pd


def main(cmdline=None):
    parser = make_parser()
    opts, args = parser.parse_args(cmdline)

    if opts.error is not None:
        return opts.error
    elif opts.bad_option:
        parser.error("you called a bad option")
    
    if len(cmdline) == 0:
        parser.print_help()
        return -1

    bam = pysam.AlignmentFile(opts.bam)
    fas = pysam.FastaFile(opts.fas)
    chr = opts.chr
    chrStart = opts.start
    chrEnd = opts.end
    outPrefix = opts.prefix

    results = pysamstats.load_variation(bam,fas,chrom=chr,start=chrStart,end=chrEnd)    
    freqA = results['A']/(results['reads_all']+1e-9)
    freqC = results['C']/(results['reads_all']+1e-9)
    freqG = results['G']/(results['reads_all']+1e-9)
    freqT = results['T']/(results['reads_all']+1e-9)
    alleleFreq = pd.DataFrame({'A':freqA,'C':freqC,'G':freqG,'T':freqT},index=results.pos)
    alleleFreq.index.name = chr
    
    alleleFreq.plot(style='.').legend(loc='center left',bbox_to_anchor=(1., 0.5))
    plt.show()


def make_parser():
    usage = """%prog args"""
    parser = OptionParser(usage)
    parser.add_option("-b","--bam",help="BAM file")
    parser.add_option("-f","--fas",help="reference genome")
    parser.add_option("-c","--chr",help="chromosome ID")
    parser.add_option("-s","--start",help="chromosome start position",type="int")
    parser.add_option("-e","--end",help="chromosome end position",type="int")
    parser.add_option("-p","--prefix",help="output prefix")
    parser.set_defaults(bad_option=False,createRDSIndex=False,error=None,make_template=False,number=0)
    return parser
 
    

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))   
