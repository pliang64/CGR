#!/usr/bin/sh
#SBATCH --account=def-pliang
#SBATCH --cpus-per-task=10
#SBATCH --time=23:59:00
#SBATCH --mem=20G
#SBATCH --nodes=1-1
#SBATCH --output=./log/%j.out

#this script requires one input as the sample ID and it generates a cram file and a bcf file and their index files
#a log folder is needed in the current folder.


if [ ! $1 ]; then echo "Usage: $0 SampleID GroupID"; exit ]; fi

module load bwa samtools samblaster bcftools

if [ ! -f $1.GCF12X.cram -a ! -s ~/rrg/lianglab-rrg/grape/CRA006917/$f.GCF12X.cram ]; then
    if [ ! -f ../$1/${1}_f1.fq.gz -o ! -f ../$1/${1}_r2.fq.gz ]; then echo "fastq files not found."; exit 1; fi
    bwa mem -t 10 -R "@RG\tID:CRA006917\tSM:$1" ~/DB/genomeSeq/Vitis_vinifera/GFC_12X ../$1/${1}_f1.fq.gz ../$1/${1}_r2.fq.gz |samblaster -r |samtools view -Sb -@ 10 -\
	|samtools sort -@ 10 -O SAM -T $1 - |samtools view -SC --reference ~/DB/genomeSeq/Vitis_vinifera/GFC_12X.fna -@ 10 -o $1.GCF12X.cram --write-index -
fi
if [ ! -s ../bcf/$1.bcf -a ! -s ~/rrg/lianglab-rrg/grape/bcf/$s.bcf ]; then
    bcftools mpileup -Ou --threads 10 -f ~/rrg/lianglab-rrg/grape/GCF_000003745.3_12X_genomic.fna $1.GCF12X.cram |bcftools call -mv -Ob --threads 10 -o ../bcf/$1.bcf
fi
if [ -s ../bcf/$1.bcf -a ! -s ../bcf/$1.bcf.csi ]; then
    bcftools index ../bcf/$1.bcf
fi

