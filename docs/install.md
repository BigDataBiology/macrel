# Install

The recommended method of installation is through
[bioconda](https://anaconda.org/bioconda/macrel):

```bash
conda install -c bioconda macrel
```

## Install from source

If you want to use an unreleased version from Github, for example, we provide a
script using _conda_ (in particular, using
[bioconda](https://bioconda.github.io/) and
[conda-forge](https://conda-forge.org/)) to install all dependencies in a
Macrel-specific environment:

```bash
git clone https://github.com/BigDataBiology/Macrel
cd macrel
./install.sh
conda activate envs/Macrel_env
```

Thereafter, to use macrel, activate this environment.

## License

Macrel is licensed under the **MIT** license (see `LICENSE` file).



