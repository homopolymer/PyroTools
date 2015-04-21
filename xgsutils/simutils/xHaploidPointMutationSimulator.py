## reference: https://www.biostars.org/p/12417/

from random import random, choice
import sys
from itertools import groupby
import datetime

def help():
    print "Simulate point mutations in given fasta sequences."
    print "Usage: xHaploidPointMutationSimulator <FASTA_FILE> <MUTATE_RATE>"
    print "-h,--help    print help message"
    exit()


def vcf_header(fasta_name):
    print "##fileformat=VCFv4.1"
    print "##fileDate=",datetime.date.today()
    print "##source=xHaploidPointMutationSimulator"
    print "##reference=",fasta_name
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTANT"


def vcf_record(record_info):
    print record_info


def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


def main(fasta_name, mutation_freq):
    for header, seq in fasta_iter(fasta_name):
        seq = list(seq)
        for i, s in enumerate(seq):
            val = random()
            if val < mutation_freq:
                ref_allele = seq[i]
                # choose a random nucleotide that's different.
                mut_allele = choice([x for x in "ACTG" if x != s.upper()])
                # vcf record
                record = "%s\t%d\t.\t%c\t%c\t.\t.\t.\tGT\t1/1"%(header,i+1,ref_allele,mut_allele)
                vcf_record(record)


if __name__ == "__main__":
    if sys.argv[1] == "-h" or sys.argv[1] == "-help" or sys.argv[1] == "--help":
        help()
    vcf_header(sys.argv[1])
    main(sys.argv[1], float(sys.argv[2]))
