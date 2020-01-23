# Install

The recommended method of installation is through
[bioconda](https://anaconda.org/bioconda/macrel):

```bash
conda install -c bioconda macrel
```

## Install from source

If you want to use an unreleased version from Github, for example, we provide a
script which _conda_ (in particular, using
[bioconda](https://bioconda.github.io/) and
[conda-forge](https://conda-forge.org/)) to install all dependencies in a
Macrel-private environment:

```bash
git clone https://github.com/BigDataBiology/Macrel
cd macrel
./install.sh
conda activate envs/Macrel_env
```

Henceforth, to use macrel, activate this environment.


## License

Macrel is licensed as GPLv3 (see `LICENSE` file).

While Macrel as a whole is **GPL v3** licensed (to comply with it being used in
some of its dependencies, namely Peptides), the macrel-specific code is also
licensed under the **MIT** license.

