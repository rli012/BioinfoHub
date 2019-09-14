#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --output=STAR_index.log
#SBATCH -p batch

#sbatch Slurm_STAR_Index_Reference.sh

#module load gatk/3.4-46
#GATK=/opt/linux/centos/7.x/x86_64/pkgs/gatk/3.4-46/GenomeAnalysisTK.jar
#Picard=/opt/linux/centos/7.x/x86_64/pkgs/picard/2.6.0/bin/picard
#JAVA=/opt/linux/centos/7.x/x86_64/pkgs/java/jdk1.7.0_17/bin/java
STAR=/rhome/rli012/bigdata/SingleCell/STAR-2.6.0a/bin/Linux_x86_64/STAR
samtools=/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools


annotation=/rhome/rli012/bigdata/SingleCell/Gencode22/anno/gencode.v22.annotation.gtf # gtf annotation file
genomeFa=/rhome/rli012/bigdata/SingleCell/Gencode22/ref/GRCh38.p2.genome.fa.noPatches # fasta sequence file
genomeDir=/rhome/rli012/bigdata/SingleCell/Gencode22/ref/ #output directory

CPU=$SLURM_NTASKS


echo "Indexing..."

$STAR --runThreadN $CPU \
	  --runMode genomeGenerate \
	  --genomeDir $genomeDir \
	  --genomeFastaFiles $genomeFa \
	  --sjdbGTFfile $annotation \
	  --sjdbOverhang 100

echo 'Done'

