This is the code running at

https://big-data-biology.org/software/macrel/

# Frontend

## Project setup
```
elm init
```

### Compiles and hot-reloads for development
```
elm reactor
```

### Compiles and minifies for production
```
elm make src/Prediction.elm --optimize
```


# Running the backend


```bash
conda create -n macrel.web python=3.8
source activate macrel.web
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge
conda install macrel uwsgi flask
```

