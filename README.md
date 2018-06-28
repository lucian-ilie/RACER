# RACER: Rapid and Accurate Correction of Errors in Reads

RACER is a C++/OpenMP program that corrects sequencing errors in high-throughput DNA sequencing data, especially designed for the Illumina platform. It does not require any additional software. 

INSTALLATION
=============================================================================

Installing with Anaconda or Bioconda:

  conda install -c mmolnar racer
  
Installing from GitHub:

  git clone https://github.com/lucian-ilie/RACER.git

Then change to the RACER directory and type:

  make

RUNNING RACER
=============================================================================

To run RACER use the command: 

  < racer > < inputReads > < correctedReads > < genomeLength > 

where 

- < racer > is the executable

- < inputReads > is the input file containin the reads; fasta or fastq

- < correctedReads > will contain the corrected reads at the end

- < genomeLength > is the approximate length of the DNA molecule that originated the reads, such as the genome length in a whole genome sequencing project 

  -if only parts of a genome were sequenced, then only the total length of those parts should be used (instead of the length of the total genome)
  
  -precise value is not necessary, an approximation will work well

###CITE
If you find RACER program useful, please cite the RACER paper:

L. Ilie, M. Molnar [RACER: Rapid and accurate correction of errors in reads](http://bioinformatics.oxfordjournals.org/content/29/19/2490) Bioinformatics, 2013
