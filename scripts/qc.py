#!/usr/bin/env python
#qc.py
#provides quality control filtering for sequencing files produced by
#next-generation sequencing platforms.
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
from MPSDemultiplexer.demultiplex import quality_control


# Fetch command-line arguments
parser = argparse.ArgumentParser(
    description='''Perform quality control filtering on reads from
    next-gen sequencing platforms''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-i','--input_fp', required=True, nargs=1,
                    type=str, help='''input filepath. Should be FASTQ format.''')
parser.add_argument('-o', '--output_fp', required=True, nargs=1, type=str,
                    help='''output filepath.''')
parser.add_argument('-O', '--out_fmt', required=True, nargs=1, type=str,
                    help='''output file format (fasta or fastq)''')
#parser.add_argument('-s', '--min_av_score', nargs=1, default=[25], type=int,
#                    help='''set the minimum average read score for a read to be
#                    accepted.''')
parser.add_argument('-t', '--min_nt_score', nargs=1, default=[28], type=int,
		    help='''set the minimum basecall score. When the average
		    score of nucleotides in a window (whose size is set with
		    the -w/--win_size parameter) drops below this value, the
		    read is truncated at the beginning of the window (i.e.,
		    subsequent nucleotides are discarded).''')
parser.add_argument('-w', '--win_size', nargs=1, default=[10], type=int,
		    help='''set the window size for determining sequence
		    truncation. A moving window scans along each read, advancing
		    one nucleotide at a time, when the average score drops below
		    that set by -t/--min_nt_score, the read is truncated.''')
parser.add_argument('-f', '--left_trim', nargs=1, default=[0], type=int,
                    help='''set the number of nucleotides to trim from
                    the front of all reads before demultiplexing.''')
parser.add_argument('-l', '--min_length', nargs=1, default=[None], type=int,
                    help='''set the minimum read length (i.e., before trimming
                    barodes, primers and any padding).''')
parser.add_argument('-L', '--max_length', nargs=1, default=[None], type=int,
                    help='''set the sequence length that should be returned
                    (i.e., after trimminmg barcodes, primers and
                    any padding).''')
#parser.add_argument('-B', '--out_buffer', nargs=1, default=[4000], type=int,
#                    help='''set the maximum number of reads to hold in memory
#                    between write operations.''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='''activate verbose mode. Indicates progress and other
                    information.''')
#Setting required variables from command line arguments.
args = parser.parse_args()
infile = args.input_fp[0]
outfile = args.output_fp[0]
outfmt = str(args.out_fmt[0])
#Setting parameters from commmand line arguments.
#average_score = args.min_av_score[0]
minimum_score = args.min_nt_score[0]
left_trimming = args.left_trim[0]
min_length = args.min_length[0]
max_length = args.max_length[0]
win_size = args.win_size[0]
verbose = args.verbose



handle1 = open(infile, "rU") # fastq sequence file
try :
    handle2 = open(outfile, "w")
except IOError:
    outdir = str(os.path.dirname(os.path.abspath(outfile)))
    os.mkdir(outdir)
    handle2 = open(outfile, "w") # create output fasta/fastq file if it
    #doesn't exist



print "\nQC run started " + strftime("%Y-%m-%d %H:%M:%S") + "."
print "Running..."

record_iter = SeqIO.parse(handle1, "fastq")

#Initializing counters:

count = 0# number of good reads
count_bad = 0# number of reads failing QC
sum_x = 0# running sum of mean read quality
sum_x_sq = 0#running sum of squared mean read quality
sum_l = 0#running sum of mean quality controlled read length
sum_l_sq = 0#running sum of squared mean quality controlled read length


#run main loop to quailty filter, assign id and trim. q_c from demultiplex
for entry in record_iter:
    record = quality_control(entry, n_trim=left_trimming,  \
                            min_score=minimum_score, \
                            min_len=min_length, win=win_size).next()
    if record:
        record = record[:max_length]
        SeqIO.write(record, handle2, outfmt)
        count += 1 #increment counter for good reads
        #x = mean(record.letter_annotations["phred_quality"])
        #l = len(record)
        #sum_x += x
        #sum_x_sq += x**2
        #sum_l += l
        #sum_l_sq += l**2
    else:
        count_bad += 1 #increment counter for bad reads
    if verbose:
        total = count + count_bad
        sys.stdout.write("Good reads: " + str(count) + ", Bad reads: " + \
                         str(count_bad))
        sys.stdout.flush()
        sys.stdout.write("\b"*len("Good reads: " + str(count) + \
                                  ", Bad reads: " + str(count_bad)))

handle1.close()
handle2.close()

if count == 0 :
    print "Error: no reads passed the quality filter"
    sys.exit(3)




#x_mean = sum_x/count #calculate the mean of mean quality score across a read
#x_std = sqrt((sum_x_sq/(count-1))-(x_mean**2)) #calculate the sample st. dev. of mean quality score
#l_mean = sum_l/count #calculate the mean length of reads
#l_std = sqrt((sum_l_sq/(count-1))-(l_mean**2)) #calculate the sample st. dev. of lengths


# Write log
logpath = open(str(os.path.splitext(outfile)[0]) + ".log","wb")
logpath.write("Logfile for quality control filtering of " \
+ infile + "\n" + strftime("%Y-%m-%d %H:%M:%S") + "\n\n" \
"Parameters specified:\n" \
"Minimum nucleotide quality score: " + str(minimum_score) + "\n" \
"Window size to check for minimum score: " + str(win_size) + "\n" \
#"Lowest acceptable average read quality score " \
#+ "(after trimming bad bases): " + str(average_score) + "\n" \
"Nucleotides to trim from front of each read: " \
+ str(left_trimming) + "\n" \
"Minimum sequence length: " + \
str(min_length) + "\n" \
"Maximum sequence length: " + \
str(max_length) + "\n\n" \
"Counts:\n" \
"Total high quality reads written: " + str(count) + "\n" \
"Total poor quality reads ignored: " + str(count_bad) + "\n")
#"Mean of mean QCed read quality: " + str(x_mean) + "\n" \
#"Standard deviation of mean QCed read quality: " + \
#str(x_std) + "\n" +\
#"Mean QCed read length: " + str(l_mean) + "\n" +\
#"Standard deviation of QCed read length: " + \
#str(l_std))
logpath.close()

print "Run finished " + strftime("%Y-%m-%d %H:%M:%S") + "."



