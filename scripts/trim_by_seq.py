#!/usr/bin/env python
#trim_by_seq.py
#trim nucleotides from 5' and 3' ends of sequences based on query strings

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

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData
import sys, os, re, argparse, resource
from time import strftime
from numpy import mean
from math import sqrt

#Add parent directory of current file to sys.path so imports can be made from folder
sys.path.insert(0,os.path.dirname(os.path.dirname(sys.argv[0])))
from MPSDemultiplexer.demultiplex import ambiguous_search, ambiguous_seq_dist


# Fetch command-line arguments
parser = argparse.ArgumentParser(
    description='''Trim nucleotides at 5' and 3' ends of a (set of) sequence(s)
    based on sequence, e.g., trim primers from front and back.''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-i','--input_fp', required=True, nargs=1,
                    type=str, help='''input filepath.''')
parser.add_argument('-o', '--output_fp', required=True, nargs=1, type=str,
                    help='''output filepath.''')
parser.add_argument('-I', '--in_fmt', required=True, nargs=1, type=str,
                    help='''input file format (fasta or fastq)''')
parser.add_argument('-O', '--out_fmt', required=True, nargs=1, type=str,
                    help='''output file format (fasta or fastq)''')
parser.add_argument('-5', '--fiveprime', required=True, nargs=1, type=str,
                    help='''5' nucleotide sequence to search for. Read(s)
                    will be truncated after the last nucleotide of this
                    sequence. Ambiguous nucleotide codes can be used.''')
parser.add_argument('-3', '--threeprime', nargs='?', type=str,
                    help='''3' nucleotide sequence to search for. Read(s)
                    will be truncated before the first nucleotide of this
                    sequence. Ambiguous nucleotide codes can be used.''')
parser.add_argument('-M', '--max_mismatch', nargs=1, default=[0], type=int,
                    help='''set the maximum acceptable number of mismathces
                    between the 5'/3' sequence and the returned match. The
                    algorithm will locate the sequence with the minimum
                    mismatch, but this parameter sets the maximum.''')
parser.add_argument('-l', '--min_len', nargs=1, default=[0], type=int,
                    help='''specify minimum sequence length after trimming.''')
parser.add_argument('-B', '--out_buffer', nargs=1, default=[4000], type=int,
                    help='''set the maximum number of reads to hold in memory
                    between write operations.''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='''activate verbose mode. Indicates progress and other
                    information.''')
#Setting variables from command line arguments.
args = parser.parse_args()
infile = args.input_fp[0]
outfile = args.output_fp[0]
infmt = str(args.in_fmt[0])
outfmt = str(args.out_fmt[0])
fiveprime = str(args.fiveprime[0])
threeprime = str(args.threeprime)
max_mismatch = args.max_mismatch[0]
min_len = args.min_len[0]
buffer_size = args.out_buffer[0]
verbose = args.verbose



handle1 = open(infile, "rU") # sequence file
try :
    handle2 = open(outfile, "w")
except IOError:
    outdir = str(os.path.dirname(os.path.abspath(outfile)))
    os.mkdir(outdir)
    handle2 = open(outfile, "w") # create output fasta/fastq file if it
    #doesn't exist



print "\nTrimming run started " + strftime("%Y-%m-%d %H:%M:%S") + "."
print "Running..."

record_iter = SeqIO.parse(handle1, infmt)

#Initializing counters:

output_buffer = []
sum_l = 0#running sum of mean quality controlled read length
sum_l_sq = 0#running sum of squared mean quality controlled read length

#run main loop to quailty filter, assign id and trim. q_c from demultiplex
for j, record in enumerate(record_iter) :
    forward_match = ambiguous_search(record, fiveprime, max_dist=max_mismatch)
    if forward_match:
        start = forward_match[1]#start
    #of new record will be end of matched 5' sequence.
    else:
        if verbose:
            print "5' sequence not found"
        start = None
    if threeprime :
        rev_match = ambiguous_search(record, threeprime, max_dist=max_mismatch)
        if rev_match:
            end = rev_match [0]#if a 3' sequence is specified, the end of the
            #new record will be the position before the start of the matched 3'
            #sequence.
        else:
            if verbose:
                print "3' sequence not found"
            end = None
    else:
        end = None
    record = record[start:end]#new record is truncated according to the
    #specified sequences.
    l = len(record)
    if l >= min_len:#test to see that trimmed sequence is above desired length
        sum_l += l
        sum_l_sq += l**2
        output_buffer.append(record)
    # write all reads if buffer is full.
    if j+1 % buffer_size == 0 :
    # write only longer list
        for i in output_buffer:
            SeqIO.write(i, handle2, outfmt)
        output_buffer = []#reset buffer
    if verbose:
        sys.stdout.write(str(j))
        sys.stdout.flush()
        sys.stdout.write("\b"*len(str(j)))

# final buffer clear        
for i in output_buffer:
    SeqIO.write(i, handle2, outfmt)
output_buffer = []

handle1.close()
handle2.close()

l_mean = sum_l/j#use final value of 'j' from main for loop to calculate mean seq length
l_std = sqrt((sum_l_sq/(j-1))-(l_mean**2))#calculate sample st. dev. of lengths.

# Write log
logpath = open(str(os.path.splitext(outfile)[0]) + ".log ","wb")
logpath.write("Logfile for sequence-based trimming of " \
+ infile + "\n" + strftime("%Y-%m-%d %H:%M:%S") + "\n\n" \
"Parameters specified:\n" \
"Input format: " + str(infmt) + "\n" \
"Output format: " + str(outfmt) + "\n" \
"5' nucleotide sequence searched for: " \
+ fiveprime + "\n" \
"3' nucleotide sequence searched for: " + \
threeprime + "\n" \
"Maximum acceptable mismatches against 5' and 3' query sequences: " + \
str(max_mismatch) + "\n\n" \
"Counts:\n" \
"Total sequences written: " + str(j) + "\n" +\
"Mean trimmed read length: " + str(l_mean) + "\n" +\
"Standard deviation of trimmed read length: " + str(l_std))
logpath.close()

print "Run finished " + strftime("%Y-%m-%d %H:%M:%S") + "."



