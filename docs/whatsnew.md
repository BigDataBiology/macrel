# What's new? (History)

## Unreleased

- Adds `query-ampsphere` command to query the AMPsphere database

## Version 1.3.0

*Released 5 December 2023*

- Adds density estimates to output
- Adds compatibility with Pyrodigal ≥3


## Version 1.2.0

*Released May 9 2022*

- Adopt pyrodigal v.0.7.3 to predict genes
- Implement a hybrid gene prediction system: training mode
  for contigs &gt; 100kbp and use pre-trained models for smaller contigs

### Bugfixes

- Fix the README files generation by the output type

## Version 1.1.0

Released November 8 2021

- Add support for bzip2 and xz compressed FASTA files
- Eliminate R dependency
- Include more extensive testing
- New feature computation code (Python implementation) made **macrel** about
  3.5 times faster than before:

| **Time \ Code** | **version 1.0.1** | **version 1.1.0** |
| :-----: | :-----: | :-----: |
| Real | 60m 44.380s | 17m 15.313s |
| User | 33m 37.604s | 9m 38.034s |
| System | 0m 16.076s | 0m 6.270s |

### Bugfixes

- Fix abundance/reads modes when single-ended reads were used

## Version 1.0.1

Released September 2 2021

### Bugfixes

- No longer crash on very short peptides (1 or 2 amino acids)

## Version 1.0.0

Released December 21 2020

### *User-visible improvements*

- Add README.md files to the output directories documenting the output files

## Version 0.6.1

Released October 29 2020

### *User-visible improvements*

- Atomicwrites is now an optional dependency. This implies a loss of
  functionality if not available (file writing is no longer atomic), but
  the atomicwrites dependency caused issues (especially on Mac)

## Version 0.6.0

Released October 10 2020

### *User-visible improvements*

- Add `--log-file` and `--log-append` command line arguments
- Add usage example in command line help message
- Output is now written atomically (_i.e._, no partial files)

## Version 0.5.0

Released May 11 2020

### User visible changes
  
- Fix bug with Prodigal by changing internal parameters. Although this is a
  bugfix, it will also change the results in some cases.
- Output table now includes the version on the header

### *Bugfixes*
  
- Fix bug with using `--force` and existing directories

## Version 0.4.0

Released March 16 2020

### *User-visible improvements*

  - Add `--keep-negatives` flag
  
  - Add warning for longer sequences

  - Remove Methionine from N-terminus
  
### *Bugfixes*
  
  - Fix for when no smORFs are present
	
## Version 0.3.1

Released January 24th 2020

### *Bugfixes*
  
- Fix for Python 2.7

## Version 0.3.0

Released January 24th 2020

### *User-visible improvements*

- Add `get-examples` command
- Add `get-smorfs` command
	  
### *Internal improvements*
  
- Convert to scikit-learn
- Updated training sets
- Fix license (must be GPL because of Peptides)

## Version 0.2.0

Released January 14th 2020

- Complete release
