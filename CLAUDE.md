# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What Macrel is

Macrel is a CLI pipeline (`macrel`) that mines antimicrobial peptides (AMPs) from
(meta)genomic data. It exposes a *subcommand interface*: `macrel COMMAND ...`, where
the command selects which stage(s) of the pipeline run. Distributed via Bioconda.

## Development commands

```bash
pip install .                 # install (models & ngl scripts ship as package_data)

python -m pytest macrel/tests # unit tests (feature computation, prediction, ampsphere)
python -m pytest macrel/tests/test_features.py::test_name   # single unit test

./run-tests.sh                # end-to-end CLI integration tests
```

`run-tests.sh` iterates over `tests/*/`, runs each `command.sh`, and diffs every
`out/macrel.out.<x>` against the checked-in `expected.<x>`. A directory named
`tests/error-*` is expected to exit non-zero. To debug/add an integration test, look at
an existing `tests/<case>/command.sh` and its `expected.*` files. **These integration
tests require the external binaries below to be on `PATH`** (installed via conda), so
they generally cannot run without a conda environment; `python -m pytest` does not need
them.

External tools (conda/bioconda), invoked via `subprocess` from `main.py`:
`ngless`, `megahit`, `paladin`, plus `pyrodigal` (imported as a Python library). MMSeqs2
and HMMER are used by the local AMPSphere query modes.

## Architecture

Everything is driven from `macrel/main.py`: `parse_args` → `validate_args` → a sequence
of `do_*` stages selected by `args.command`. The commands compose the same stages:

- **peptides** (amino-acid FASTA) → `do_predict`
- **contigs** (nucleotide FASTA) → `do_smorfs` → `do_predict` → `do_density`
- **reads** (FASTQ) → `do_assembly` (trim + megahit) → `do_smorfs` → `do_predict` → `do_density`
- **abundance** (FASTQ + peptide ref) → `do_read_trimming` → paladin map → ngless count
- **get-smorfs** (nucleotide FASTA) → `do_smorfs` only
- **query-ampsphere** → `do_ampsphere_query` (see below)
- **get-examples** → downloads example data, returns early (no tempdir)

The prediction core (`do_predict`) is where the science lives and is fully in-Python:

1. `AMP_features.fasta_features` computes a 22-column feature table per peptide.
   The actual physicochemical calculations are in `macrel_features.py` (adapted from
   modlAMP), with amino-acid scales/constants in `database.py`. `normalize_seq` strips a
   leading `M` and trailing `*` before feature computation.
2. `AMP_predict.predict` loads two gzipped **ONNX** classifiers via `onnxruntime`
   (`data/models/AMP.onnx.gz` = is-AMP, `Hemo.onnx.gz` = hemolytic), thresholds
   probabilities at 0.5, and derives the AMP family label (ADP/ALP/CDP/CLP) from
   charge and cysteine content. Non-AMPs are dropped unless `--keep_negatives`.

smORF prediction (`ORFs_prediction.predict_genes`) uses **pyrodigal**; the code guards
for both the pre-3.0 `OrfFinder` and the newer `GeneFinder` class name.
`filter_smorfs.filter_smorfs` keeps only peptides ≤100 aa and (with `--cluster`)
dedupes identical sequences.

Stages run inside a single `tempfile.TemporaryDirectory`; intermediate files
(trimmed reads, SAM, all-ORFs FASTA) live there. Final outputs go to `args.output` as
gzipped, tab-separated tables prefixed with a `# Prediction from macrel v<version>`
comment line. Each command also writes a `README.md` into the output dir (templates in
`output.py`). Output writing uses `utils.open_output` (atomic writes).

### AMPSphere querying

`ampsphere.py` queries the AMPSphere database in three modes (`--query-mode
exact|mmseqs|hmmer`). Default is the remote HTTP API (throttled with `sleep`); `--local`
downloads the database once (`database.py` handles caching under `--cache-dir`) and runs
matching locally with MMSeqs2/HMMER.

## Conventions

- **Version** is single-sourced in `macrel/macrel_version.py`; `pyproject.toml`
  reads it via `[tool.setuptools.dynamic]` (`version = {attr = ...}`).
  Bump it there and add a `ChangeLog` entry when releasing.
- Slow imports (`onnxruntime`, `AMP_features`, pandas-heavy modules) are imported
  *inside* the `do_*` functions, not at module top level, to keep CLI startup fast.
  Preserve this pattern.
- The `.ngl` scripts in `data/scripts/` are NGLess pipelines (trimming, counting)
  shipped as package data and invoked as external `ngless` scripts.
- CI (`.github/workflows/python-app.yml`) runs both `pytest` and `run-tests.sh` across
- Commits should go to the `main` branch; releases are tagged with `v<version>`
  and pushed to PyPI and Bioconda.

  Linux/macOS/Windows and Python 3.10–3.13.
