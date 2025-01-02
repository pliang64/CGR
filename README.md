# CGR
A Data Processing Pipeline for Whole-genome Sequencing-Based Grapevine Genetic Testing using Chaos Game Representation and Integrative Deep Learning Approaches

Here we provide a list of scripts and python codes for processing grapevine whole genome sequence (WGS) data to be used grapevine species and cultivar identification using a combination of Chaos Game Representation (CGR) and Integrative Deep Learning Approaches. The untitlity of individual script/code is described below:

1. fastq2cram2vcf4graham.sbatch: a slurm script for aligning WGS data in Illumina PE reads to a reference genome and generate cram and bcf files using high performacing computing.
2. VCF2FaSeq_v3.pl: a perl script for generating a pseudo fasta sequence to represent the genotype of variants of WGS samples for use in sequence similarity and phylogeny ananlysis
3. A set of python codes for converting fasta sequences into CGR images and for four predictive models.
