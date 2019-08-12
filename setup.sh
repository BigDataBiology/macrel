#!/bin/bash

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

read -p "[ Installing message ] :: Do you have already installed in your path conda? (1 - yes; 0 - no)
>>>>>> " conda_if

if [[ $conda_if == "1" ]]
then
	echo "[ ... Skipping initial procedures ]"
elif [[ -n $conda_if ]]
then
	echo "[ Getting Bioconda -- It can take a while... ]"
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
	echo "[ User ignored the question, sorry. Closing! ]"
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
install.packages("Peptides", dependencies=TRUE)
##########################################################################
" > inst.R 
chmod +x inst.R
R --vanilla --slave < inst.R
rm -rf inst.R

echo "# Closing environment"
conda deactivate

echo "[ ## 2.] Installing prodigal_modified_version"

echo "# Downloading from private repository"
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/bitmap.c?download=1' --output 'bitmap.c'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/bitmap.o?download=1' --output 'bitmap.o'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/bitmap.h?download=1' --output 'bitmap.h'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/training.o?download=1' --output 'training.o'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/training.h?download=1' --output 'training.h'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/training.c?download=1' --output 'training.c'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/sequence.o?download=1' --output 'sequence.o'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/sequence.h?download=1' --output 'sequence.h'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/sequence.c?download=1' --output 'sequence.c'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/node.o?download=1' --output 'node.o'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/node.h?download=1' --output 'node.h'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/node.c?download=1' --output 'node.c'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/metagenomic.o?download=1' --output 'metagenomic.o'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/metagenomic.c?download=1' --output 'metagenomic.c'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/metagenomic.h?download=1' --output 'metagenomic.h'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/dprog.o?download=1' --output 'dprog.o'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/dprog.c?download=1' --output 'dprog.c'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/dprog.h?download=1' --output 'dprog.h'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/main.h?download=1' --output 'main.h'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/main.c?download=1' --output 'main.c'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/main.o?download=1' --output 'main.o'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/gene.o?download=1' --output 'gene.o'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/gene.h?download=1' --output 'gene.h'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/gene.c?download=1' --output 'gene.c'
curl --header 'Host: zenodo.org' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:68.0) Gecko/20100101 Firefox/68.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://zenodo.org/record/3365790' --cookie 'session=2ef0bf5258eb1d00_5d51015b.2g1V17pWfFh_vy4zSPcgGsfN6-U; __atuvc=1%7C29%2C0%7C30%2C0%7C31%2C8%7C32%2C3%7C33; __atssc=google%3B2; _pk_id.57.a333=88b4a4e2ffa2de44.1562733750.12.1565591273.1565589267.; _pk_ref.57.a333=%5B%22%22%2C%22%22%2C1565589267%2C%22https%3A%2F%2Fmail.google.com%2Fmail%2Fu%2F0%2F%22%5D; __atuvs=5d50ff12d27a3709002; __atrfs=ab/|pos/|tot/|rsi/5d48059d00000000|cfc/|hash/0|rsiq/|fuid/fc115b49|rxi/|rsc/|gen/1|csi/|dr/; _pk_ses.57.a333=*' --header 'Upgrade-Insecure-Requests: 1' 'https://zenodo.org/record/3365790/files/Makefile?download=1' --output 'Makefile'
echo "# Unpacking"
make install
rm -rf *.o *.h *.c Makefile
chmod +x prodigal
mv prodigal envs/FACS_env/bin/

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
