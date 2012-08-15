#!/usr/bin/env python
#separate_assigned_reads.py
#splits a FASTA file containing reads from multiple samples into spearate files
#where each file contains reads from a single sample. It does so based on the
#names in the fasta headers for each read.
#    Copyright (C) <2012>  <Benjamin C. Smith>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#To run the program, enter "separate_assigned_reads.py inputfile mappingfile"
#on the command line. Output files will be written to a folder called
#"seperate" in the directory containing the input file.


            
from Bio import SeqIO
from collections import deque
import sys, os, argparse, resource

#resource.setrlimit(resource.RLIMIT_NOFILE, (2048,-1))

parser = argparse.ArgumentParser(
    description='''Split a fasta file containing reads from multiple samples
    by sample.''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-i','--input_fp', required=True, nargs=1,
                    type=str, help='''input filepath.''')
parser.add_argument('-I', '--in_fmt', required=True, nargs=1, type=str,
                    help='''input file format (fasta or fastq)''')
parser.add_argument('-O', '--out_fmt', required=True, nargs=1, type=str,
                    help='''required output file format (fasta or fastq)''')
parser.add_argument('-m', '--map_fp', required=True, nargs=1, type=str,
                    help='''mapping filepath. Mapping file must be a text
                    file in TSV format with a column of sample names as the
                    first column.''')
parser.add_argument('-t', '--trim_header', action='store_true', 
                    help='''enables trimming of additional info from header
                    line leaving just the piece before the first space.''')
parser.add_argument('-v', '--verbose', action='store_true', 
                    help='''enables verbose mode. Records will be printed to
                    the screen as they are written.''')
parser.add_argument('-a', '--aligned', action='store_true', 
                    help='''specifies that input file should be alignment.
                    If enabled, sequence lengths will be checked to ensure
                    they are all the same.''')
args = parser.parse_args()
infile = args.input_fp[0]
mapfile = args.map_fp[0]
trim = args.trim_header
infmt = args.in_fmt[0]
outfmt = args.out_fmt[0]
verbose = args.verbose
aligned = args.aligned

handle1 = open(infile, "rU") # input file. fasta format
handle2 = open(mapfile, "rU") # mapping file, tsv format
if aligned:
    outpath = str(os.path.dirname(os.path.abspath(infile))) \
    + "/separated_alignments/"
else:
    outpath = str(os.path.dirname(os.path.abspath(infile))) \
    + "/separated/"
try:
    os.mkdir(outpath)
except OSError:
    pass
    
record_iter = SeqIO.parse(handle1,infmt)

outfiles = {}
for line in handle2 :
    if line[0] != "#" :
        separated = line.split()
        outfiles[separated[0]] = open(outpath + str(separated[0]) + \
        ".fasta", 'wb')

lengths = deque()
for i, rec in enumerate(record_iter) :
    #Progress indicator
    if verbose:
        sys.stdout.write("Number Separated: " + str(i))
        sys.stdout.flush()
        sys.stdout.write("\b"*len("Number Separated: " + str(i)))
    lengths.append(len(rec.seq))
    if trim :
        rec.description = ""
    if aligned and min(lengths) != max(lengths) :
        print "Error raised by sequence " + rec.id
        print "Alignments do not all have the same length."
        print "Either run without -a flag or check alignment file."
        sys.exit(1)
    SeqIO.write(rec, outfiles[rec.id.split("_")[0]], outfmt)

for keys in outfiles:
    outfiles[keys].close()
handle1.close()
handle2.close()

