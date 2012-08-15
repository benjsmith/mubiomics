#!/usr/bin/env python
#seqparser.py
#Provides various sequence manipulation tools.
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
import sys, argparse

parser = argparse.ArgumentParser(
    description='''Can produce the reverse, reverse compliment or compliment of
    sequences in a multiple sequence file. Can also convert between any file
    formats supported by Biopython's SeqIO class (e.g., FASTQ to FASTA and
    FASTA + QUAL to FASTQ).''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-i','--input_fp', required=True, nargs=1,
                    type=str, help='''input filepath. Should be a FASTA format file.''')
parser.add_argument('-o', '--out_fp', required=True, nargs=1, type=str,
                    help='''output filepath. .''')
parser.add_argument('-q', '--qual_fp', required=False, nargs='?', type=str,
                    help='''(optional) input qual filepath for converting
                    FASTA + QUAL to FASTQ.''')
parser.add_argument('-r', '--rev_comp', action='store_true', 
                    help='''produces the reverse compliment of the input
                    sequences.''')
parser.add_argument('-c', '--comp', action='store_true', 
                    help='''produces the compliment of the input sequences.''')
parser.add_argument('-e', '--rev', action='store_true', 
                    help='''produces the reverse of the input sequences.''')
parser.add_argument('-I', '--in_fmt', required=True, nargs=1, type=str,
                    default='fasta',
                    help='''specify the format of input sequences.''')
parser.add_argument('-O', '--out_fmt', required=True, nargs=1, type=str,
                    default='fasta',
                    help='''specify the format of output sequences.''')

args = parser.parse_args()
infile = args.input_fp[0]
outfile = args.out_fp[0]
qualfile = args.qual_fp
rev_comp = args.rev_comp
comp = args.comp
rev = args.rev
infmt=args.in_fmt[0]
outfmt=args.out_fmt[0]

handle1 = open(infile, "rU")
try :
    handle2 = open(outfile, "wb")
except IOError:
    outdir = str(os.path.dirname(os.path.abspath(outfile)))
    os.mkdir(outdir)
    handle2 = open(outfile, "wb") # output fasta file


if infmt=='fasta' and outfmt=='fastq':
    try:
        qualhandle = open(qualfile, "rU")
        record_iter = SeqIO.QualityIO.PairedFastaQualIterator(handle1, qualhandle)
    except TypeError:
        print "Cannot convert from FASTA to FASTQ without QUAL file. \
See \"seq_parser.py --help\" for details."
        sys.exit(1)
else:
    record_iter = SeqIO.parse(handle1,infmt)


if comp and rev_comp:
    print "Error: Cannot set -c(--comp) and -r(--rev_comp) in same run."
    sys.exit(1)

print "Conversion running..."

for rec in record_iter:
    new_seq = rec
    if comp :
        new_seq.seq = Seq(str(rec.seq.complement()), IUPAC.unambiguous_dna)    
    if rev_comp :
        new_seq.seq = Seq(str(rec.seq.reverse_complement()), IUPAC.unambiguous_dna)
    if rev :
        new_seq.seq = Seq(str(rec.seq)[::-1], IUPAC.unambiguous_dna)
    SeqIO.write(new_seq, handle2, outfmt)


handle1.close()
handle2.close()
try:
    qualhandle.close()
except NameError:
    pass
print "Conversion finished."
