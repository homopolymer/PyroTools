
## Introduction

**PyroTools** is a toolkit including methods developed for the analysis of human genome sequencing data and metagenomic sequencing data.

##Methods 

- **StrainBinning**

  A method for binning long reads or assembled contigs into the strain-level groups, and resolving the issue of strain heterogeneity in the coassembly of multiple strains.  Our preliminary results show that it could assemble the complete or most complete individual genomes of closely related strains from metagenomic sequencing data, when using the isolated data binned by StrainBinning.  Comparatively, existing assemblers commonly assemble sequencing data of closely related strains into short and chimeric contigs.  StrainBinning internally uses StrainCall.
  
  [1] Zeng F, Chen T. (2015) Strain-level binning of long reads for metagenomics.  Unpublished. 

- **StrainCall** 

  A method for reconstructing local strain sequences from metagenomic sequencing data.  It is ultra-fast, and accurate in both the estimation of compositional abundances and the reconstruction of strain sequences. It works for both shotgun and amplicon sequencing data.  StrainCall accurately explores the intraspecies diversity in a local region, and provides a foundation for purifying metagenomics data that mixs up multiple strains and allowing the isolation of long reads/assembled contigs into the strain-level groups.  It eventually facilitates the complete assembly of individual strains, as what we have done in StrainBinning.

  [2] Zeng F, Chen T. (2015) StrainCall: a fast method for reconstructing strain sequences from metagenomic sequencing data. Submitted.

- **ProbAlign** 
  
  A program to estimate the parameters of the alignment scoring function, and to correct the misalignments.
  
  [3] Zeng F, Jiang R, Ji G, Chen T. (2014) ProbAlign: a re-alignment method for long sequencing reads. [_biorXiv.org_](http://biorxiv.org/content/early/2014/09/02/008698).
  
- **SnpCall/PyroHMMsnp**
 
  A program to call SNPs in the sequencing data of an individual.

  [4] Zeng F, Jiang R, Chen T. (2013) PyroHMMsnp: a SNP caller for 454 and Ion Torrent sequencing data. _Nucleic Acids Research_.

- **IndelCall/PyroHMMvar**

  A program to call insertions and deletions in the sequencing data of an individual.

  [5] Zeng F, Jiang R, Chen T. (2013) PyroHMMvar: a sensitive and accurate method to detect short Indels and SNPs for Ion Torrent and 454 data. _Bioinformatics_.



## Clone and compile
	> git clone https://github.com/homopolymer/PyroTools.git
	> cd PyroTools
	> chmod 775 run_install.sh
	> ./run_install.sh 
	
	
### Dependencies
Require the following packages:

- [Samtools](samtools.sourceforge.net)
- [Bamtools](https://github.com/pezmaster31/bamtools)
- [Seqtk](https://github.com/lh3/seqtk)
- [Bedtools](bedtools.readthedocs.org)
- [NGSUtils](ngsutils.org)
- [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) 
- [Bowtie2](bowtie-bio.sourceforge.net/bowtie2)








