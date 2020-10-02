# What's new? (History)

### *User-visible improvements*

- Add `--log-file` flag

### *User-visible improvements*

- Add usage example in command line help message

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
