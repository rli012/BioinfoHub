### STAR: https://github.com/alexdobin/STAR
wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz	
tar -xzf 2.7.3a.tar.gz  
cd STAR-2.7.3a  
cd source  
make STAR  
	
### samtools: http://www.htslib.org/download/
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2  
tar -xjf samtools-1.9.tar.bz2  
cd samtools-1.9  
./configure --prefix=/home/rli3/bin/samtools-1.9  
make  
make install  
	
**vim ~/.bashrc   
export PATH=/home/rli3/bin/STAR-2.7.3a/bin/Linux_x86_64:$PATH  
export PATH=/home/rli3/bin/samtools-1.9/bin:$PATH**  
	
### Subread: http://subread.sourceforge.net/
Download a Subread binary distribution that suits your oprating system  
wget https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz  
tar xvzf subread-2.0.0-Linux-x86_64.tar.gz  
	
### FastQC: https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip  
unzip fastqc_v0.11.8.zip  
cd FastQC  
chmod 755 fastqc  
	
### fastp: https://github.com/OpenGene/fastp
this binary was compiled on CentOS, and tested on CentOS/Ubuntu  
wget http://opengene.org/fastp/fastp  
chmod a+x ./fastp  
OR  
get source (you can also use browser to download from master or releases)  
git clone https://github.com/OpenGene/fastp.git  
cd fastp  
make  
	
### MultiQC: https://github.com/ewels/MultiQC
virtualenv env  
source env/bin/activate  
export PYTHONPATH=`pwd`/env/lib/python2.7/site-packages  
\#pip install --upgrade pip  
pip install multiqc  
	
**scp -r FastQC/ rli012@cluster.hpcc.ucr.edu:~/bigdata/G/**  
