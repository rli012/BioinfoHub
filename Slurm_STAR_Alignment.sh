#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --output=alignment.log
#SBATCH -p batch


#module load gatk/3.4-46
STAR=/rhome/rli012/bigdata/SingleCell/STAR-2.6.0a/bin/Linux_x86_64/STAR
#GATK=/opt/linux/centos/7.x/x86_64/pkgs/gatk/3.4-46/GenomeAnalysisTK.jar
#Picard=/opt/linux/centos/7.x/x86_64/pkgs/picard/2.6.0/bin/picard
#JAVA=/opt/linux/centos/7.x/x86_64/pkgs/java/jdk1.7.0_17/bin/java
samtools=/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools

annotation=/rhome/rli012/bigdata/SingleCell/anno/gencode.v19.annotation.gtf
genomeDir=/rhome/rli012/bigdata/SingleCell/ref/
genomeFa=/rhome/rli012/bigdata/SingleCell/ref/GRCh37.p13.genome.fa

N=$SLURM_ARRAY_TASK_ID
CPU=$SLURM_NTASKS

FILE=`ls raw/SRR*\.fastq.gz | grep _1.fastq.gz | head -n $N | tail -n 1`
PREFIX=${FILE%_1.fastq.gz}
PREFIX=${PREFIX#raw/}

fq1=$FILE
fq2=${FILE/_1/_2}

#PREFIX=SRR2973290
#bam=${PREFIX}.bam
#fq1=${PREFIX}_1.fastq.gz
#fq2=${PREFIX}_2.fastq.gz

echo 'Start Alignment...'
echo $PREFIX

### Alignment ###

$STAR --runThreadN $CPU \
	  --genomeDir $genomeDir \
      --twopassMode Basic \
	  --readFilesIn $fq1 $fq2 \
	  --readFilesCommand zcat \
	  --outSAMtype BAM Unsorted \
      --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
      --clip5pNbases 5 \
      --clip3pNbases 5 \
      --limitBAMsortRAM 19732153018 \
	  --outFileNamePrefix alignment/${PREFIX}

$samtools sort -@ $CPU alignment/${PREFIX}Aligned.out.bam -T alignment/${PREFIX} -o alignment/${PREFIX}.bam
rm alignment/${PREFIX}Aligned.out.bam


echo 'Done!'

