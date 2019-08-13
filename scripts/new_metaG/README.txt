README

REQUIREMENTS

Perl (with standard packages: GetOpt, File and Math installed)
BioPerl (with SeqIO and Seq installed)

if you cannot used precompiled executables:
GCC
as well as standard bash utils ar and ranlib

argtable2 library (is part of the sourcecode zip file)

INSTALLATION

1. download files

2a. if you need to compile the executables:
#extract sourcecode
tar xvfz simulatorCparts.tar.gz 
#change directory
cd simulatorCparts/ 
#compile
make clean
make 

2a. copy the executables to a path of your choice and set this path as "$binfolder" in the the main simulator perl file.

3. run program
perldoc simulatorForSolexaReads.pl
	for infos and you will see following:

==> Different version of the program:
==> simulatorForSolexaReads_1.2.pl: updated version of the original
==> simulatorForSolexaReads_1.2_exact.pl: version without stochastic abundance estimation; abundances are directly transformed into readcounts
==> simulatorForSolexaReads_1.2_exact_cp.pl: same as simulatorForSolexaReads_1.2_exact.pl but using convert_project from the mira assembler (needs to be edited for the location of your local convert project_copy, but is much faster then BioPerl)


NAME
       iMESSi - Simulator for Illumina metagenomic sequences
SYNOPSIS
       simulatorForSolexaReads.pl [options]
	   
OPTIONS
--help					Prints this help message
--outfolder				Folder for output
--outprefix				Prefix of output filenames (default: SolexaSimReads)
--abundanceFile			input file with genomes and their organism number
--genomeInfo			file with information about the genomes: genome length, copy number in an organism.
--genomeFolder			folder the genomes mentioned in abundanceFile are at.
--readlength			Read length of simulated reads.
--insertNumber			Number of inserts to simulate, from each insert 2 paired-end reads are generated.
--insertSize			Mean insert size of simulated reads.
--insertSD				Standard deviation from the mean insert size.

--qualityfiles			location of quality files for sequencing error and quality generation. 
Note that the first file can be given by just specifying it's location and multiple files can be given using folling scheme --qualityfiles="file1.qual -q file2.qual -q file3.qual". [Sample files can be found at http://www.bork.embl.de/~mende/simulated_data/*.qual.gz]

--binfolder				Location of generateReads and simulateSequencingErrors, has to be set here or inside the script

			
Options to add more data to old simulations
--multiplexIdentifier   Identifier for multiplex Libraries (default is 0)
--laneNumber            Starting lane number
--tileNumber            Starting tile number
--xCoord                Starting xCoord number
--yCoord                Starting yCoord number

DESCRIPTION
This program is used to simulate Illumina reads from a specified genome or metagenome
	   
AUTHOR
Daniel Mende <mende@embl.de>, Alison Waller,  et al. 
Please cite us: 
Mende DR, Waller AS, Sunagawa S, Järvelin AI, Chan MM, et al. (2012) Assessment of Metagenomic Assembly Using Simulated Next Generation Sequencing Data. PLoS ONE 7(2): e31386. doi:10.1371/journal.pone.0031386

We like to thank Falk Hildebrand for implementing the automatic genome info file generation part and Paul Igor Costea for improving the C-code
This program is open source, feel free to send in improved versions.

VERSIONS
Version 1.0 - First release
Version 1.1 - Fixed some bugs and added automatic genome info file generation
Version 1.2 - Added file that gives the exact location every read originates from
Version 1.3 - Improved C-code with high speed-up.
