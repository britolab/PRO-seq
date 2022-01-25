
# installation notes

# source $HOME/miniconda3/bin/activate
# export PATH=/home/acv46/.local/bin:$PATH
# wget https://raw.githubusercontent.com/LaboratorioBioinformatica/MARVEL/master/environment.yml
# conda env create -n marvel -f=environment.yml

source $HOME/miniconda3/bin/activate
conda activate marvel
export OMP_NUM_THREADS=8

marvel --help

conda deactivate
