#!/bin/bash
#$ -N Download_ArrayExpress
# #$ -V # Job has the same environment variables as the submission shell
#$ -l h=xxx08
#$ -pe peomics 10 # -pe threads
# #$ -l mem_free=8G # h_vmem=4G hard limit of the maximum amount of vitual memory
#$ -o "logs/Download_ArrayExpress.stdout"
#$ -e "logs/Download_ArrayExpress.stderr"
#$ -wd "/home/xxx"
# #$ -cwd

# #$ -v "REF=hg38,FQ=test.fq,SAM=test.sam" # no space!!!
# #$ -l h_rt=00:02:00
# #$ -wd "/home/rli3/Documents"
# #$ -t 1-3 # task array

# usage: qsub SGE_qsub_template.sh
# OR
# usage: qsub -N run_test -V -pe threads 2 -l h_vmem=4G -o "logs/test.stdout" -e "logs/test.stderr" SGE_qsub_template.sh

CPU=$NSLOTS

echo 'Start'

FILE=`cut -f34 /home/data/ArrayExpress/E-MTAB-xxxx/E-MTAB-xxxx.sdrf.txt | grep 'ftp' | tail -n +8`

for fl in $FILE; do
	echo $fl
	`wget $fl`;
done

echo 'Done!'

