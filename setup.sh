#!/usr/bin/env bash

echo "
#########################################################################
######## FACS pipeline - Fast AMP Clustering System   ###################
#########################################################################
########                    Authors:                  ###################
########                                              ###################
######## Célio Dias Santos Júnior, Luis Pedro Coelho  ###################
#########################################################################
######## ISTBI - FUDAN University / Shanghai - China  ###################
#########################################################################

"

if ! which conda > /dev/null; then
	echo "[ Conda not found. Getting Bioconda -- It can take a while... ]"
	curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	
	sh Miniconda3-latest-Linux-x86_64.sh
	
	export PATH=$PATH:~/miniconda3/bin/
	eval "$(conda shell.bash hook)"
	
	conda config --add channels r	
	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge

	conda config --set auto_activate_base false
else
	echo "[ Conda found. ]"
fi

export PATH=$PATH:~/miniconda3/bin/
eval "$(conda shell.bash hook)"
echo "[ ## 1.] Installing routine linux softwares"
echo "# Creating environment of execution of routine softwares"
mkdir envs
conda create -p ./envs/FACS_env python=3.7
conda activate ./envs/FACS_env

echo "# Installing routine linux softwares"
conda install -y sqlite
conda install -y trimmomatic
conda install -y megahit
conda install -y pandaseq
conda install -y paladin
conda install -y samtools
conda install -y eXpress
conda install -y pigz
conda install -y parallel

echo "# Installing routine python packages"
conda install -y matplotlib
conda install -y scikit-learn
conda install -y pandas

echo "# Installing routine R packages"
conda install -y r-essentials
conda install -y r-base
conda install -y r-caret
conda install -y r-randomforest
conda install -y r-dplyr
conda install -y r-data.table

echo "#!/usr/bin env

##########################################################################
install.packages(\"Peptides\", repos = \"http://cran.us.r-project.org\", dependencies=TRUE)
install.packages(\"doParallel\", repos = \"http://cran.us.r-project.org\", dependencies=TRUE)
install.packages(\"obliqueRF\", repos = \"http://cran.us.r-project.org\", dependencies=TRUE)
##########################################################################
" > inst.R 
chmod +x inst.R
R --vanilla --slave < inst.R
rm -rf inst.R

echo "# Closing environment"
conda deactivate

echo "[ ## 2.] Installing prodigal_modified"

(cd prodigal_modified && make && mv prodigal ../envs/FACS_env/bin/prodigal_sm )

echo "[ ## 3.] Getting python scripts"
echo "# Downloading"
curl -O https://raw.githubusercontent.com/Superzchen/iLearn/master/descproteins/CTDDClass.py
curl -O https://raw.githubusercontent.com/Superzchen/iLearn/master/descproteins/readFasta.py 
curl -O https://raw.githubusercontent.com/Superzchen/iLearn/master/descproteins/saveCode.py
echo "# Modifying"
chmod +x *.py

echo "############ Installation procedures finished
****** Thank you for installing FACS ********
--- If any bug appears report to:
celio.diasjunior@gmail.com

"
