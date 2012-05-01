#!/usr/bin/env python
#fasta2fastq.py
#Fasta/qual to fastq converter. 
 
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

# Fetch command-line arguments
parser = argparse.ArgumentParser(
    description='''Convert paired fasta and qual files to a single fastq
    file.''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-f','--fasta_fp', required=True, nargs=1,
                    type=str, help='''input fasta filepath.''')
parser.add_argument('-q', '--qual_fp', required=True, nargs=1, type=str,
                    help='''input qual filepath.''')
parser.add_argument('-o', '--output_fp', required=True, nargs=1, type=str,
                    help='''output fastq filepath.''')
parser.add_argument('-B', '--buffer_size', nargs=1, default=[4000], type=int,
                    help='''set the maximum number of reads to hold in memory
                    between write operations.''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='''activate verbose mode. Indicates progress and other
                    information.''')



#Setting required variables from command line arguments.
args = parser.parse_args()
fastafile = args.fasta_fp[0]
qualfile = args.qual_fp[0]
outfile = args.output_fp[0]
buffer_size = args.buffer_size[0]
verbose = args.verbose

print strftime("%Y-%m-%d %H:%M:%S") + " - Running..."

fastahandle = open(fastafile, "rU")
qualhandle = open(qualfile, "rU")
record_iter = SeqIO.QualityIO.PairedFastaQualIterator(fastahandle, qualhandle)
output_buffer = []

try :
    outhandle = open(outfile, "wb")
except IOError:
    outdir = str(os.path.dirname(os.path.abspath(outfile)))
    os.mkdir(outdir)
    outhandle = open(outfile, "wb") # output fasta file

output_buffer = []
for i, record in enumerate(record_iter):
    if verbose:
        sys.stdout.write(str(i))
        sys.stdout.flush()
        sys.stdout.write("\b"*len(str(i)))       
    output_buffer.append(record)
    if i % buffer_size == 0:
        for seq in output_buffer:
            SeqIO.write(seq, outhandle, "fastq")
        output_buffer = []

#final buffer flush
for seq in output_buffer:
    SeqIO.write(seq, outhandle, "fastq")
output_buffer = []

print strftime("%Y-%m-%d %H:%M:%S") + " - Finished."

fastahandle.close()
qualhandle.close()
outhandle.close()
