# WEBSERVER

This is the code running at

http://big-data-biology.org/software/macrel/


# RUNNING THE WEBSERVER


```bash
conda create -n macrel.web python=3.7
source activate macrel.web
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge
conda install macrel uwsgi flask
```
