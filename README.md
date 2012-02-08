# Massively-Parallel Sequencing Software

Authors: Benjamin C. Smith and David M. Moskowitz, 2011

## Description:

A package for processing short-read sequences generated by
second-generation sequencing platforms (e.g. Roche 454, Illumina HiSeq2000).

Code is pre-publication and in development, but mostly works. Please cite us if you use it
([email](mailto:benjamin.smith@einstein.yu.edu) for citation information)

## Contents:

- `demultiplexer.py` a sequence demultiplexer for assigning reads from
multiple samples based on DNA barcode, written in python. Currently
supports 8-bp Hamming barcodes and regular barcodes of any fixed
length.

- `qc.py` a quality control script. Filters and trims reads based on a
number of settable parameters.

- `trim_by_seq.py` a script to trim, e.g., primers from both ends of
reads in a sequence file. Allows specification of a 5' sequence
and an optional 3' reverse complement sequence. Accepts ambiguous
nucleotide values in the 5' and 3' sequences. Returned reads will begin
after the 5' sequence and end before the 3' sequence. If either sequence
is not located, that end will not be trimmed.

- `ucstripper.py`  a python script for processing the output from a
usearch --cluster run against a reference database. The output is an OTU
table where each column corresponds to a sample and each row corresponds
to sequence in the database matched by a read. The cells are filled with
count data.

- `ucstripper_paired.py`  identical to `ucstripper.py` except that it's
designed to work with reads that are from a paired-end run. Reads can
be analysed in usearch as though they are separate, but this script will
count a hit to an OTU only if both ends agree on the assigned OTU. If
only one end is present in the .uc file, its OTU is counted. 

- `pool_otus.py` a python script which takes an OTU table and pools
counts whose OTU names are identical. Any OTUs that didn't find a match
in the databse are combined into a category called "Noise". It also
allows setting of a minimum count threshold; if no sample contains more
counts tnan this threshold, the OTU is discarded.

For all included python scripts, type the name of the script followed by
`-h` or `--help` to see help documentation.

- shell scripts for splitting input and output to distribute processes over
multiple cores (currently only suitable for the qc.py script).

## Installation:

To install, place the mubiomics directory anywhere on your hard drive, add it to the $PYTHONPATH shell variable and add the mubiomics/scripts directory to your $PATH variable.

The following prerequisites must be installed for scripts to work:
-Python v2.7
-Biopython v1.5.8 or greater
-matplotlib v1.0, or greater, and its dependencies (e.g. numpy, scipy etc.).

Some of the scripts are designed to process results from the following read
classification programs. To run the shel scripts that use them, they should
also be installed:
-[usearch](http://www.drive5.com/usearch/)
-[rdp multiclassifier](http://rdp.cme.msu.edu/classifier/classifier.jsp)
-[pplacer](http://matsen.fhcrc.org/pplacer/)

## Testing:
To run tests, cd into mpsdemultiplexer/tests and enter the following on the command line,
followed by return:

    $ tests.sh
	
In order for the tests to run to completion usearch v4.0, or greater must be installed.
Once downloaded, the binary must be placed on the shell $PATH.

If the test fails, you will get error messages which may indicate the reason for the failure(s).
It's most likely to be either that the dependencies are not installed and/or configured correctly
or that the programs are not on your shell $PATH variable. To see $PATH, type `echo $PATH`
in terminal. To add programs to your shell $PATH and $PYTHONPATH, enter the following
in `~/.bash_profile`:

    export PATH=${PATH}:/path/to/directory
    export PYTHONPATH=${PYTHONPATH}:/path/to/directory

The `tests.sh` script also provides an example workflow. Open it in a text editor for explanations of each step.


## QIIME compatibility:

The programs in this package were intially written so that they provided
compatibility with [QIIME](http://qiime.sourceforge.net/), which, at the time, couldn't demultiplex FASTQ data from the Hi-Seq platform. Output from
demultiplexer.py is thus compatible with all downstream QIIME workflows.
Output from ucstripper.py and pool\_otus.py are an identical format to
the OTU tables produced by QIIME and thus compatible with the analysis scripts that take them as input.

## Notes:
Development and testing was performed on Mac OSX 10.6 using Python v2.7, numpy v1.5.1, maplotlib v1.0.1, Biopython v1.5.8. We can't guarantee that it will work with other setups, but feel free to [email](mailto:benjamin.smith@einstein.yu.edu) with any issues.

The patricia.py class was obtained from a [post on stack overflow](http://stackoverflow.com/questions/2406416/implementing-a-patricia-trie-for-use-as-a-dictionary/2412468#2412468). Thank you to John Peel for posting this, we hope you don't mind us using it!