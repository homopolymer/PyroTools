PyroTools
===============

**Sincere apology to users** that PyroHMMsnp and PyroHMMvar have been out of maintenance for a period because the important change of my life and job.  God blesses.  Everything is fine at the end.  And then, we are on the way to integrate these two tools into the toolkit, named PyroTools.  The code is changed frequently and will be launched soon.

## Introduction

**PyroTools** is a toolkit for the processing of *next-generation sequencing* data, e.g. adjust the inconsistent alignments in BAM file(s), and analysis of next-generation sequencing data, e.g. detect SNPs and short InDels.

**Methods** in PyroTools include

- ProbAlign, a program to estimate the parameters of the alignment scoring function, and to correct the misalignments.
- SnpCall, a program to call SNPs in the sequencing data of an individual \[1\].
- IndelCall, a program to call insertions and deletions in the sequencing data of an individual \[2\].


**Features** of PyroTools include

-  BAM misalignment correction 
-  Sequencing error modelling
-  Haplotype based variant calling
-  Feasible computational time


## Re-alignment

### A example

- The original mappings of 454 WGS data of NA12878 \[3\]

![Resize icon][ex2_rawdata]
[ex2_rawdata]:./realignment_examples/ex2_rawdata.png

- The re-alignment results of PyroTools

![Resize icon][ex2_pyrotools]
[ex2_pyrotools]:./realignment_examples/ex2_PyroTools.png

- The re-alignment results of [SRMA](sourceforge.net/projects/srma/)

![Resize icon][ex2_srma]
[ex2_srma]:./realignment_examples/ex2_SRMA.png

### More examples

- Example 1.  The top track is the original mappings.  The bottom track is the re-alignment results of PyroTools.

![Resize icon][more_ex1]
[more_ex1]:./realignment_examples/ex4_chr21_11022500_11022530.png

- Example 2. 

![Resize icon][more_ex2]
[more_ex2]:./realignment_examples/ex6_chr21_11032590_11032620.png

- Example 3.

![Resize icon][more_ex3]
[more_ex3]:./realignment_examples/ex12_chr21_11052670-11052700.png

- Example 4.

![Resize icon][more_ex4]
[more_ex4]:./realignment_examples/ex13_chr21_11061190-11061240.png


===================================================================

\[1\]: Zeng F., Jiang R., Chen T. PyroHMMsnp: a SNP caller for 454 and Ion Torrent sequencing data. _Nucleic Acids Research_, 2013.

\[2\]: Zeng F., Jiang R., Chen T. PyroHMMvar: a sensitive and accurate method to detect short Indels and SNPs for Ion Torrent and 454 data. _Bioinformatics_, 2013.

\[3\]: The BAM files of 454 WGS data of NA12878 are available in 1KG project ftp site.





