#!/bin/bash
#$ -N mm10
#$ -V # Job has the same environment variables as the submission shell
#$ -l h=xxx12
#$ -pe peomics 10 # -pe threads
# #$ -l mem_free=8G # h_vmem=4G hard limit of the maximum amount of vitual memory
#$ -o "logs/FastQC.stdout"
#$ -e "logs/FastQC.stderr"
#$ -wd "/home/rli3/xxx/data/"
# #$ -cwd

# #$ -v "REF=hg38,FQ=test.fq,SAM=test.sam" # no space!!!
# #$ -l h_rt=00:02:00
# #$ -wd "/home/rli3/Documents"
# #$ -t 1-3 # task array

# usage: qsub SGE_qsub_template.sh
# OR
# usage: qsub -N run_test -V -pe threads 2 -l h_vmem=4G -o "logs/test.stdout" -e "logs/test.stderr" SGE_qsub_template.sh

fastqc=~/bin/FastQC/fastqc
CPU=$NSLOTS
FILES=`ls fastq/ERR*\.fastq.gz`

$fastqc $FILES --outdir=FastQC/

