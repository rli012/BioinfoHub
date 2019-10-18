#!/bin/bash
#$ -N mm10
#$ -V # Job has the same environment variables as the submission shell
#$ -l h=xxx16
#$ -pe peomics 10 # -pe threads
# #$ -l mem_free=8G # h_vmem=4G hard limit of the maximum amount of vitual memory
#$ -o "logs/fastp.QC.stdout"
#$ -e "logs/fastp.QC.stderr"
#$ -wd "/home/rli3/Projects/xxx/data"
# #$ -cwd

# #$ -v "REF=hg38,FQ=test.fq,SAM=test.sam" # no space!!!
# #$ -l h_rt=00:02:00
# #$ -wd "/home/rli3/Documents"
# #$ -t 1-3 # task array

# usage: qsub SGE_qsub_template.sh
# OR
# usage: qsub -N run_test -V -pe threads 2 -l h_vmem=4G -o "logs/test.stdout" -e "logs/test.stderr" SGE_qsub_template.sh


CPU=$NSLOTS

FILES=`ls fastq/ERR*\.fastq.gz`

for FILE in $FILES
do
	PREFIX=${FILE%.fastq.gz}
	PREFIX=${PREFIX#fastq/}

	fq1=$FILE
	fqo1=${fq1/fastq/trimmed_fastq}
	fqo1=${fqo1/fastq./trimmed.fastq.}
#fq2=${FILE/_1/_2}

#PREFIX=SRR2973290
#bam=${PREFIX}.bam
#fq1=${PREFIX}_1.fastq.gz
#fq2=${PREFIX}_2.fastq.gz

	echo 'Start Alignment...'
	echo $PREFIX
	echo $fqo1
### fastp QC ###

	fastp -i $fq1 -o $fqo1 -h fastpQC/${PREFIX}.fastp.html -j fastpQC/${PREFIX}.fastp.json
done

echo 'Done!'


