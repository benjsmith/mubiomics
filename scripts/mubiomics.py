#!/usr/bin/env python

#Paired-end demultiplexing program

#Author: Benjamin Smith, 2011

from Bio import SeqIO
import sys, argparse, os
import cPickle as pickle
from time import strftime

sys.path.insert(0,os.path.dirname(os.path.dirname(sys.argv[0])))
from MPSDemultiplexer.demultiplex import *

parser = argparse.ArgumentParser(
    description='''A set of python tools for characterizing the microbiome
    from next-generation sequencing reads.''',
    fromfile_prefix_chars='@')
parser.add_argument('-i','--infile', required=True, nargs=1,
                    type=str, help='''input filepath. Should be a FASTA format
                    file. Usually will be the file of unassigned sequences
                    produced by the demultiplexer.''')
parser.add_argument('-I', '--in_fmt', required=True, nargs=1, type=str,
                    help='''input file format (fasta or fastq)''')
parser.add_argument('-t','--matefile', required=True, nargs=1,
                    type=str, help='''filepath of mate file to be demultiplexed.
                    Should be a FASTA or FASTQ format file. The mates of all
                    in the input file must be in here and in the same order,
                    though they do not neccessarily have to appear
                    sequentially in the mate file. This will be the default
                    situation if you have used the demultiplexer in this
                    package.''')
parser.add_argument('-T', '--mate_fmt', required=True, nargs=1, type=str,
                    help='''mate file format (fasta or fastq)''')
parser.add_argument('-d', '--direction_ids', required=True, nargs=2,
                    type=int, help='''numbers identifying direction of mate.
                    Two numbers, separated by a space, should be entered,
                    e.g., if forward reads end in "/1" and reverse reads end
                    in "/3", you should enter -d 1 3''')
parser.add_argument('-o', '--outfile', required=True, nargs=1, type=str,
                    help='''output filepath. Demultiplexed mate file will be
                    written here.''')
parser.add_argument('-O', '--out_fmt', required=True, nargs=1, type=str,
                    help='''required output file format (fasta or fastq)''')
parser.add_argument('-m', '--map_fp', nargs=1, type=str, default=[None],
                    help='''(optional) mapping filepath.
                    Should be TSV format containing a list of sequences
                    to be searched for that should be trimmed from the front
                    of the mate sequences before being written, e.g.,
                    reverse primers. The primer name should be in the first
                    column and the nucleotide sequence in the seccond column.
                    Comment lines in this file should begin with a "#"''')
parser.add_argument('-M', '--max_mismatch', nargs=1, default=[3], type=int,
                    help='''set the maximum acceptable number of mismathces
                    for the sequences listed in the mapping file
                    that will be trimmed. If no mapping file is sepcified
                    this value will be ignored.''')
parser.add_argument('-f', '--left_trim', nargs=1, default=[None], type=int,
                    help='''(optional) set the number of nucleotides to trim from the
                    front of the mate sequence before writing. This can be
                    used to trim primers if all are of the same length.''')
parser.add_argument('-L', '--max_length', nargs=1, default=[None], type=int,
                    help='''set the sequence length that should be returned
                    (i.e., after any trimminmg that is specified).''')
parser.add_argument('-d', '--use_indexdb', action='store_true',
                    help='''activate database indexing mode. This should be
                    used if the mate file is larger than the available RAM. It
                    indexes the entries in the mate file and temporarily
                    stores them in an SQL database file.''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='''activate verbose mode to indicate progress.''')
#set parameters from supplied arguments
args = parser.parse_args()
infile = args.infile[0]
matefile = args.matefile[0]
matefmt = args.mate_fmt[0]
infmt = args.in_fmt[0]
outfmt = args.out_fmt[0]
dir_ids = args.direction_ids
outfile = args.outfile[0]
mapfile = args.map_fp[0]
max_mismatch = args.max_mismatch[0]
trim = args.left_trim[0]
max_length = args.max_length[0]
use_indexdb = args.use_indexdb
verbose = args.verbose

#open file handles
#inhandle = open(infile, 'rU')
#matehandle = open(matefile, 'rU')
if mapfile:
    maphandle = open(mapfile, "rU")
    trimlist = [str(line.split()[1]) for line in maphandle if line[0]!="#"]
#get data
#produce dictionary-like objects where sequences in the files are indexed by
#their id
print "\nMate finding run started " + strftime("%Y-%m-%d %H:%M:%S") + "."
indata = SeqIO.parse(open(infile, "rU"), infmt)
#open output file for writing
try :
    outhandle = open(outfile, "wb")
except IOError:
    outdir = str(os.path.dirname(os.path.abspath(outfile)))
    os.mkdir(outdir)
    outhandle = open(outfile, "wb")
    
print "Indexing mate sequences..."

if use_indexdb:
    indexfile = str(os.path.splitext(os.path.abspath(matefile))[0]) + ".idx"
    matedata = SeqIO.index_db(indexfile, matefile, matefmt)
else:
    matedata = SeqIO.index_db(matefile, matefmt)