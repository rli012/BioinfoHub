
### Server Environment
module unload miniconda2
module load miniconda3
module unload python2
module load python/3.6.0


### Initialize (Run once)
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge # do not add to channels?

conda init bash && source ~/.bashrc


### Configure (in ~/.condarc)
channels:
  - defaults
  - bioconda
pkgs_dirs:
  - ~/bigdata/.conda/pkgs
envs_dirs:
  - ~/bigdata/.conda/envs
auto_activate_base: false


### Creat Conda Environment
conda info
conda create -n envconda python=3.6.7


### Activating
conda activate envconda


### Installing Packages
conda install -n envconda whatshap


### Deactivating
conda deactivate


### SBATCH
active the virtual envrioment by running *conda activate envconda*, then submit the job as usual

