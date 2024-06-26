# Install

The recommended method of installation is through
[bioconda](https://anaconda.org/bioconda/macrel):

```bash
conda install -c bioconda macrel
```

## Install from source

If you want to use an unreleased version from Github, or tweak the code, we still recommend conda to manage the dependencies. You can install from source by running:

```bash
git clone https://github.com/BigDataBiology/Macrel
cd macrel

conda create -n macrel_env
conda activate macrel_env
conda install -y \
          ngless \
          pyrodigal \
          megahit \
          paladin \
          pandas \
          "scikit-learn<1.3.0" \
          "joblib<1.3.0" \
          atomicwrites \
          tzlocal

pip install .
```

Thereafter, to use macrel, activate this environment.

## License

Macrel is licensed under the **MIT** license (see `LICENSE` file).



