#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --output=featureCounts.log
#SBATCH -p batch


#module load gatk/3.4-46

featureCounts=/rhome/rli012/bigdata/software/subread-2.0.0-Linux-x86_64/bin/featureCounts

annotation=/rhome/rli012/bigdata/SingleCell/anno/gencode.v19.annotation.gtf
genomeDir=/rhome/rli012/bigdata/SingleCell/ref/
genomeFa=/rhome/rli012/bigdata/SingleCell/ref/GRCh37.p13.genome.fa

N=$SLURM_ARRAY_TASK_ID
CPU=$SLURM_NTASKS

BAMS=`ls alignment/*\.bam`
OUTPUT='featureCounts.txt'

echo 'Start Counting...'

### Alignment ###

$featureCounts -T $CPU --primary --ignoreDup -g gene_name -a $annotation -o count.tmp $BAMS

tail -n+2 count.tmp | cut --complement -f2-5 > $OUTPUT

echo 'Done!'
