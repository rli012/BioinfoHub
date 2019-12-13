#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=24G
#SBATCH --time=5:00:00
#SBATCH --output=Check_number_of_reads.log
#SBATCH -p batch
#SBATCH --chdir=/bigdata/jialab/rli012/PCa/data/fromSRA/GSE54460/


#sbatch --array 1-101 download_sra.pair.sh

module load sratoolkit/2.10.0


CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=$1
fi

FILES=`cat SRR_Acc_List.txt`

for FILE in $FILES
do
	echo $FILE
	FQ1=$FILE\_1.fastq.gz
	FQ2=$FILE\_2.fastq.gz

	zcat $FQ1 | wc -l
	zcat $FQ2 | wc -l
done
