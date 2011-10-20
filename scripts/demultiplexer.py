#!/usr/bin/env python

#Demultiplexing program for Illumina reads. 
 
#Author: Benjamin Smith, 2011

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData
import sys, os, re, argparse, resource
from time import strftime
from numpy import *

#Add parent directory of current file to sys.path so imports can be made from folder
sys.path.insert(0,os.path.dirname(os.path.dirname(sys.argv[0])))
from MPSDemultiplexer.demultiplex import *
from MPSDemultiplexer.patricia import *
from MPSDemultiplexer.hamming import *



# Fetch command-line arguments
parser = argparse.ArgumentParser(
    description='''Demultiplex barcoded reads from
    next-gen sequencing platforms, trimming primers and barcodes
    and producing FASTA format files containing assigned and unassigned
    sequences, respectively, and a log file. Optionally can handle paired
    end reads where the barcode is only on one end. The mate read names
    with no barcode will be prefixed with the same name as the barcoded
    read and will appear at the same position in the mated output file.''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-i','--input_fp', required=True, nargs=1,
                    type=str, help='''input filepath.''')
parser.add_argument('-o', '--output_fp', required=True, nargs=1, type=str,
                    help='''output filepath.''')
parser.add_argument('-p', '--pair_fp', nargs='?', type=str,
                    help='''(optional) filepath to mate file of a paired-end
                    sequencing run. If specified, output will consist of an
                    additional FASTA file, suffixed with ".mates". The
                    ordering of demultiplexed reads in the two output
                    sequence files will be the same. All reads that had no
                    mates will be placed in the FASTA file suffixed with
                    ".unassigned".''')
parser.add_argument('-I', '--in_fmt', required=True, nargs=1, type=str,
                    help='''input file format (fasta or fastq)''')
parser.add_argument('-O', '--out_fmt', required=True, nargs=1, type=str,
                    help='''required output file format (fasta or fastq)''')
parser.add_argument('-m', '--map_fp', required=True, nargs=1, type=str,
                    help='''mapping filepath.''')
parser.add_argument('-L', '--max_length', nargs=1, default=[None], type=int,
                    help='''set the sequence length that should be returned
                    (i.e., after trimminmg barcodes, primers and
                    any padding).''')
parser.add_argument('-b', '--barcode_type', nargs=1, default=['hamming_8'], type=str,
                    help='''set the barcode type. Supported types are:
                    'hamming_8', or a number speifying the barcode length.
                    If a number is given, specifying max_barcode_errors > 0
                    (-e option) may be risky, as sequencing errors could
                    lead to false matching. Using Hamming codes reduces the
                    chances of this significantly, by correcting single
                    nucleotide errors.''')
parser.add_argument('-M', '--max_primer_mismatch', nargs=1, default=[3], type=int,
                    help='''set the maximum acceptable number of mismathces
                    in the primer sequences. If no primers in sequences, e.g.,
                    if they are fragmented DNA, set this number equal to
                    the max. acceptable barcode errors.''')
parser.add_argument('-e', '--max_barcode_errors', nargs=1, default=[0], type=int,
                    help='''set the maximum acceptable number of barcode errors.
                    ''')
parser.add_argument('-T', '--max_bc_start_pos', nargs=1, default=[0], type=int,
                    help='''set the maximum number of nucleotides in the read
                    before the barcode would be expected to start.
                    E.g., if there are 3 bp of padding to the left of the
                    barcode, this should be set to 3.''')
parser.add_argument('-R', '--right_pad', nargs=1, default=[0], type=int,
                    help='''set the number of nucleotides separating
                    the end of the barcode and the beginning of
                    the target primer sequence.''')
parser.add_argument('-n', '--start_numbering_at', nargs=1, default=[1], type=int,
                    help='''set the start number for labeling assigned reads''')
parser.add_argument('-S', '--separate', action='store_true', 
                    help='''use this if a separate FASTA file is required
                    for each sample (barcode). These files will be
                    written to a new directory called "separated", created
                    within the directory of the output filepath. A single
                    file containing all assigned sequences will still be
                    created at the output filepath even if this option
                    is specified. This can be deleted by the user if it is
                    not required.''')
parser.add_argument('-B', '--out_buffer', nargs=1, default=[4000], type=int,
                    help='''set the maximum number of reads to hold in memory
                    between write operations.''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='''activate verbose mode. Indicates progress and other
                    information.''')
#Setting required variables from command line arguments.
args = parser.parse_args()
infile = args.input_fp[0]
infmt = args.in_fmt[0]
pairfile = args.pair_fp
outfmt = args.out_fmt[0]
outfile = args.output_fp[0]

# tab-delimited meta-data file
mapfile = args.map_fp[0]

#Setting parameters from commmand line arguments.
max_length = args.max_length[0]
barcode_type = args.barcode_type[0]
try :
    barcode_length = int(barcode_type)
except ValueError :
    if barcode_type == 'hamming_8' :
        barcode_length = 8
    else :
        print "--barcode_type (-b) specified is an unsupported type"
        sys.exit(1)
max_bc_mismatch = args.max_barcode_errors[0]
max_prim_mismatch = args.max_primer_mismatch[0] 
max_bc_st_pos = args.max_bc_start_pos[0]
right_padding = args.right_pad[0]
out_buffer = args.out_buffer[0]
verbose = args.verbose


inhandle = open(infile, "rU") # fastq sequence file
try :
    outhandle = open(outfile, "wb")
    if pairfile:
        mateouthandle = open(str(os.path.splitext(outfile)[0]) + \
".mates.fasta", "wb") # mate output fasta file
except IOError:
    outdir = str(os.path.dirname(os.path.abspath(outfile)))
    os.mkdir(outdir)
    outhandle = open(outfile, "wb") # output fasta file
    if pairfile:
        mateouthandle = open(str(os.path.splitext(outfile)[0]) + \
".mates.fasta", "wb") # mate output fasta file
maphandle = open(mapfile, "rU") # mapping file, tsv format
uahandle = open(str(os.path.splitext(outfile)[0]) + \
".unassigned.fasta", "wb") # unassigned output fasta file


print "\nDemultiplexing run started " + strftime("%Y-%m-%d %H:%M:%S") + "."

record_iter = SeqIO.parse(inhandle, infmt)


bcs = {} #Initialize bcs dictionary. Will have barcode as keyword 
#and position in list as value
ids = [] #Initialize list of sample IDs
primers = [] #Initialize list of primers
line_ct = -1
#Extract barcodes, primers and IDs from mapping file
for line in maphandle :
    if line[0] != "#" :
        line_ct += 1
        separated = line.split()
        ids.append(separated[0])
        bcs[separated[1]] = line_ct
        try:
            primers.append(separated[2])
        except IndexError:
            pass
#Test that mapping file yielded some data.
try :
    ids[0]
except IndexError:
    print "No data loaded from mapping file. Check and try again."
    sys.exit(1)

print "Building radix trie " + strftime("%Y-%m-%d %H:%M:%S") + "..."
bc_trie = mmDNAtrie(bcs, max_bc_mismatch) # from demultiplex
print "Built " + strftime("%Y-%m-%d %H:%M:%S") + "."
print "Running..."

#Initializing counters:
# id_counts = matrix storing number of reads in each ID from mapping file [0]
# count = number of reads
# count_bad = number of reads failing QC
# count_u = unassigned reads
# count_a = assigned reads [initial value for output formatting of number of reads in category]
# count_pm = sequences with primer mismatches, including unassigned that were rejected for too many primer mismatches
# primer_mismatches = number of primer mismatches in each read
# primer_lengths = length of primer of each read
# point_error = primer_mismatches/primer_lengths
# ua_output_buffer = unassigned reads output buffer
# a_output_buffer = assigned reads output buffer

id_counts = zeros( (len(ids), 2), dtype=object)
for i, id in enumerate(ids) :
    id_counts[i, 0] = id

count = 0
count_bad = 0
count_u = 0
count_a = args.start_numbering_at[0]-1
count_pm = 0
primer_mismatches = [] 
primer_lengths = [] 
point_error = [] 
ua_output_buffer = [] 
a_output_buffer = []

#run main loop to assign id and trim.
for record in record_iter :
    result = identify_read(bc_trie, record, bcs, ids, primers, \
            max_pos=max_bc_st_pos, max_p_mismatch=max_prim_mismatch, \
            rpad=right_padding, bc_len=barcode_length, bc_type=barcode_type) # from demultiplex
    if max_length:
        max_trim_pos = result[0]+max_length
    else:
        max_trim_pos = None
    if result[1] != "Unassigned" :
        new_seq_r = record[result[0]:max_trim_pos]
        new_seq_r.description = "%s orig_bc=%s new_bc=%s bc_diffs=%i" \
                        % (new_seq_r.name, result[2], result[3], \
                        ambiguous_seq_dist(result[2], result[3]))                    
        id_counts[bcs[result[3]],1] += 1
        count_a += 1
        primer_mismatches.append(float(result[4]))
        primer_lengths.append(float(result[5]))
        if result[5] != 0:
            point_error.append(float(result[4])/float(result[5]))
        else:
            point_error.append(0)
        new_seq_r.id = "%s_%i" % (result[1], count_a)
#             print new_seq_r.format("fasta").rstrip("\n\b")
        a_output_buffer.append(new_seq_r)
    else:
        count_u += 1
        new_seq_r = record[result[0]:max_trim_pos]
        ua_output_buffer.append(new_seq_r)
    if result[4]:
        count_pm += 1
    
    
    # write all for assigned or unassigned reads if buffer is full.
    if count_u + count_a >= out_buffer :
        # write only longer list
        if len(ua_output_buffer) > len(a_output_buffer):
            for i in ua_output_buffer:
                SeqIO.write(i, uahandle, outfmt)
            ua_output_buffer = [] #reset unassigned buffer
        else:
            for i in a_output_buffer:
                SeqIO.write(i, handle2, outfmt)
            a_output_buffer = [] #reset assigned buffer
    count += 1 #increment counter for reads    
#Progress indicator
    if verbose:
        sys.stdout.write("Assigned: " + str(count_a) + ", Unassigned: " + str(count_u))
        sys.stdout.flush()
        sys.stdout.write("\b"*len("Assigned: " + str(count_a) + ", Unassigned: " + str(count_u)))
# final buffer clear        
for i in ua_output_buffer:
    SeqIO.write(i, uahandle, outfmt)
ua_output_buffer = []
for i in a_output_buffer:
    SeqIO.write(i, outhandle, outfmt)
a_output_buffer = []

if len(primer_mismatches) > 0:
    av_pm = mean(primer_mismatches)
    sd_pm = std(primer_mismatches)
    tot_pm = sum(primer_mismatches)
    av_err = mean(point_error)
    sd_err = std(point_error)
else:
    av_pm = 0
    sd_pm = 0
    tot_pm = 0
    av_err = 0
    sd_err = 0

logpath = open(str(os.path.splitext(outfile)[0]) + ".log ","wb")
logpath.write("Logfile for demultiplexing of " \
+ infile + "\n" + strftime("%Y-%m-%d %H:%M:%S") + "\n\n" \
"Parameters specified:\n" \
"Maximum sequence length: " + \
str(max_length) + "\n" \
"Maximum barcode mismatches: " + str(max_bc_mismatch) + "\n" \
"Maximum primer mismatches: " + str(max_prim_mismatch) + "\n" \
"Maximum position in read at which barcode could start: " + \
str(max_bc_st_pos) + "\n" \
"Number of bases between end of barcode and start of primer " \
+ str(right_padding) + "\n\n" \
"Counts:\n" \
"Total reads parsed: " + str(count) + "\n" \
"Assigned reads written: " + str(count_a-args.start_numbering_at[0]+1) + "\n" \
"Reads exceeding acceptable mismatches: " + str(count_pm) + "\n" \
"Unassigned reads: " + str(count_u) + "\n" \
"(Discard stats. below if no primers were specified) \n" \
"Average primer errors in assigned reads: " + str(round(av_pm, 4)) + "\n" \
"Std. dev. of primer errors in assigned reads: " + str(round(sd_pm, 4)) + "\n" \
"Total primer errors for all assigned reads: " + str(int(tot_pm)) + "\n" \
"Estimated error rate in assigned reads: " + str(round(av_err, 4)) + "\n" \
"Std. dev. of estimated error rate in assigned reads: " + \
str(round(sd_err, 4)) + "\n\n" \
"Sample id.\tCounts\n")
savetxt(logpath, id_counts, fmt='%s', delimiter='\t')

inhandle.close()
outhandle.close()
maphandle.close()
uahandle.close()
logpath.close()

#If user specified separate output files for each sample (barcode), the
#following code is invoked.
if args.separate :
    print "Separating...\n"
    sephandle = open(outfile, "rU")
    maphandle = open(mapfile, "rU")
    try:
        seppath = str(os.path.dirname(os.path.abspath(outfile))) \
        + "/separated/"
    except OSError:
        print "No input file specified, two expected."
        sys.exit(1)
    try:
        os.mkdir(seppath)
    except OSError:
        pass

    sep_iter = SeqIO.parse(sephandle,"fasta")
    resource.setrlimit(resource.RLIMIT_NOFILE, (1500,-1))
    outfiles = {}
    for line in maphandle :
        if line[0] != "#" :
            ls = line.split()
            outfiles[ls[0]] = open(seppath + str(ls[0]) + \
            ".fasta", 'wb')

    for i, rec in enumerate(sep_iter) :
        if verbose:
            sys.stdout.write("Number Separated: " + str(i))
            sys.stdout.flush()
            sys.stdout.write("\b"*len("Number Separated: " + str(i)))
#        if i < 1000 :
        SeqIO.write(rec, outfiles[rec.id.split("_")[0]], "fasta")
#        else:
#            break
    
    for keys in outfiles:
        outfiles[keys].close()
    sephandle.close()
    maphandle.close()


print "Run finished " + strftime("%Y-%m-%d %H:%M:%S") + "."



