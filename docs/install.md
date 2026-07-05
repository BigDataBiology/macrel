# Install

The recommended method of installation is through
[bioconda](https://anaconda.org/bioconda/macrel), using
[pixi](https://pixi.sh/):

```bash
pixi global install -c conda-forge -c bioconda macrel
```

You can also install it with conda (or mamba):

```bash
conda install -c bioconda macrel
```

## Run without installing

If you have [pixi](https://pixi.sh/) installed, you can run Macrel directly
without installing it first:

```bash
pixi exec -c conda-forge -c bioconda macrel
```

This downloads Macrel (and its dependencies) into a temporary environment and
runs it, which is convenient for one-off usage.

## Install from source

If you want to use an unreleased version from Github, or tweak the code, we still recommend conda to manage the dependencies. You can install from source by running:

```bash
git clone https://github.com/BigDataBiology/Macrel
cd macrel

conda create -n macrel_env
conda activate macrel_env
conda install -y -c conda-forge -c bioconda \
          ngless \
          pyrodigal \
          megahit \
          paladin \
          pandas \
          requests \
          onnxruntime \
          atomicwrites

pip install .
```

The `-c conda-forge -c bioconda` channels are required: `ngless`, `megahit`,
`paladin`, and `pyrodigal` are distributed on bioconda and will generally not
resolve from the default channels.

Thereafter, to use macrel, activate this environment.

## Optional: local AMPSphere querying

Querying the AMPSphere database locally (`macrel query-ampsphere --local`) with
the `mmseqs` or `hmmer` modes requires MMSeqs2 and/or HMMER to be installed:

```bash
conda install -c bioconda mmseqs2 hmmer
```

## License

Macrel is licensed under the **MIT** license (see `LICENSE` file).



