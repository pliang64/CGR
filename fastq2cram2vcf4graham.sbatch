#!/usr/bin/sh
#SBATCH --account=[yourUserID]
#SBATCH --cpus-per-task=10
#SBATCH --time=23:59:00
#SBATCH --mem=20G
#SBATCH --nodes=1-1
#SBATCH --output=./log/%j.out

#this script requires one input as the sample ID and it generates a cram file and a bcf file and their index files
#required folders at the current location: "log", "crams", "bcf"
#the input fastq.gz sequence files for each sample, $1, need to be named as "$1.1.fastq.gz" and "$1.2.fastq.gz" in the ./fastq folder
$the output files including the corresponding cram files named as "$1.GCF12X.cram" and "$1.GCF12X.cram.crai" in ./cram
#as well as bcf files in ./bcf folder as "$1.bcf" and "$1.bcf.csi"
#group ID $2 is optional

if [ ! $1 ]; then echo "Usage: $0 SampleID GroupID"; exit ]; fi
if  [ ! $2 ]; then $2=$1; fi

module load bwa samtools samblaster bcftools

if [ ! -f ./crams/$1.GCF12X.cram ]; then
    if [ ! -f ./fastq/${1}.1.fastq.gz -o ! -f ./fastq/${1}.1.fastq.gz  ]; then echo "fastq files not found."; exit 1; fi
    bwa mem -t 10 -R "@RG\tID:$2\tSM:$1" GFC_12X /fastq/${1}.1.fastq.gz /fastq/${1}.2.fastq.gz |samblaster -r |samtools view -Sb -@ 10 -\
	|samtools sort -@ 10 -O SAM -T $1 - |samtools view -SC --reference GCF_12X.fna -@ 10 -o crams/$1.GCF12X.cram --write-index -
fi
if [ ! -s ../bcf/$1.bcf -a ! -s ./bcf/$s.bcf ]; then
    bcftools mpileup -Ou --threads 10 -f GCF_12X.fna ./crams/$1.GCF12X.cram |bcftools call -mv -Ob --threads 10 -o ./bcf/$1.bcf
fi
if [ -s ./bcf/$1.bcf -a ! -s ./bcf/$1.bcf.csi ]; then
    bcftools index ./bcf/$1.bcf
fi
