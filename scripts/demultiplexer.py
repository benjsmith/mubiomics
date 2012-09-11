#!/usr/bin/env python
#demultiplexer.py
#Demultiplexing program for next-gen sequencing reads. 
 
#
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
#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
#from Bio.Data import IUPACData
import sys, os, re, argparse, resource
from time import strftime
#from numpy import *
#from collections import deque, Counter
#from pympler.summary import summarize

#Add parent directory of current file to sys.path so imports can be made from folder
sys.path.insert(0,os.path.dirname(os.path.dirname(sys.argv[0])))
from MPSDemultiplexer.demultiplex import *
from MPSDemultiplexer.patricia import *
from MPSDemultiplexer.hamming import *



# Fetch command-line arguments
parser = argparse.ArgumentParser(
    description='''Demultiplex barcoded reads from
    next-gen sequencing platforms, trimming primers and barcodes
    and producing sequence files containing assigned and unassigned
    sequences, respectively, and a log file. Optionally, it can handle paired
    end reads where the barcode is only on one end. The mated read names
    from the end with no barcode will be prefixed with the same name as the
    barcoded end read and will appear at the same position in the mated output
    file. WARNING: This can be very slow or use a large amount of
    memory, depending on which options you specify and on the sizes of your input
    sequence files. If barcodes are on both ends, run each file
    separately, as though demultiplexing for single end reads. If barcodes on
    each ends are different, you could also join the two sequence files and
    demultiplex it in one, using the mapping file to control naming of forward
    and reverse reads.
    ''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-i','--input_fp', required=True, nargs=1,
                    type=str, help='''input filepath.''')
parser.add_argument('-o', '--output_fp', required=True, nargs=1, type=str,
                    help='''output filepath.''')
parser.add_argument('-p', '--pair_fp', nargs='?', type=str,
                    help='''(optional) filepath to mate file of a paired-end
                    sequencing run. If specified, output will consist of an
                    additional file, suffixed with ".mates". The
                    ordering of demultiplexed reads in the two output
                    sequence files will be the same. All demultiplexed reads
                    that had no mates will be placed in the output file suffixed
                    with ".singletons". All reads that could not be identified
                    to a sample name will be written to the output file suffixed
                    with ".unassigned".''')
parser.add_argument('-P', '--paired', action='store_true',
                    help='''activate paired-end mode.''')
parser.add_argument('-H', '--hardware', required=True, nargs=1,
                    help='''enter the sequencing hardware used. Options are
                    '454', 'HiSeq', 'MiSeq'.''')
parser.add_argument('-I', '--in_fmt', required=True, nargs=1, type=str,
                    help='''input file format (fasta or fastq)''')
parser.add_argument('-O', '--out_fmt', required=True, nargs=1, type=str,
                    help='''required output file format (fasta or fastq)''')
parser.add_argument('-m', '--map_fp', required=True, nargs=1, type=str,
                    help='''mapping filepath. This should be a tab-separated
                    file. See example in mubiomics/tests/testmap.txt. If
                    demultiplexing paired ends, the reverse primer should be
                    in column 4, i.e., the next column after the forward
                    primer.''')
parser.add_argument('-L', '--max_length', nargs=1, default=[None], type=int,
                    help='''set the sequence length that should be returned
                    (i.e., after trimminmg barcodes, primers and
                    any padding).''')
parser.add_argument('-b', '--barcode_type', nargs=1, default=['hamming_8'],
                    type=str, help='''set the barcode type. Supported types are:
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
parser.add_argument('-s', '--suppress_unassigned', action='store_true',
                    help='''activate this mode to stop unassigned sequences
                    being written. This can save space if you don't care about
                    retaining sequences that couldn't be assigned an I.D.''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='''activate verbose mode. Indicates progress and other
                    information.''')
#Setting required variables from command line arguments.
args = parser.parse_args()
infile = args.input_fp[0]
infmt = args.in_fmt[0]
pairfile = args.pair_fp
paired = args.paired
#dir_ids = args.direction_ids
hardware=args.hardware[0]
outfmt = args.out_fmt[0]
outfile = args.output_fp[0]
max_length = args.max_length[0]
barcode_type = args.barcode_type[0]
# tab-delimited meta-data file
mapfile = args.map_fp[0]

try :
    barcode_length = int(barcode_type)
except ValueError :
    if barcode_type == 'hamming_8' :
        barcode_length = 8
    else :
        print "--barcode_type (-b) " + str(barcode_type) + " is an unsupported type"
        print "Run demultiplexer.py -h for help"
        sys.exit(1)
max_bc_mismatch = args.max_barcode_errors[0]
max_prim_mismatch = args.max_primer_mismatch[0] 
max_bc_st_pos = args.max_bc_start_pos[0]
right_padding = args.right_pad[0]
#use_indexdb = args.use_indexdb
#idx_exists = args.index_exists
suppress = args.suppress_unassigned
verbose = args.verbose

# Set behaviour for determinig read ids and read direction from
# hardware and dir_ids parameters.
if hardware == "454":
    split_on = ''
    rn_part = 0
elif hardware == "HiSeq":
    split_on = '/'
    rn_part = 0
elif hardware == "MiSeq":
    split_on = ' '
    rn_part = 0
else:
    print "Unknown sequencing hardware specified (-H, --hardware, parameter)."
    print "Run demultiplexer.py -h for help"
    sys.exit(1)


# Open required filehandles dependent on parameters specified

try :
    outhandle = open(outfile, "wb")
    if paired:
        mateouthandle = open(str(os.path.splitext(outfile)[0]) + \
".mates." + outfmt, "wb") # mate output fasta file
        singouthandle = open(str(os.path.splitext(outfile)[0]) + \
".singletons." + outfmt, "wb") # singleton output fasta file
except IOError:
    outdir = str(os.path.dirname(os.path.abspath(outfile)))
    os.mkdir(outdir)
    outhandle = open(outfile, "wb") # output fasta file
    if paired:
        mateouthandle = open(str(os.path.splitext(outfile)[0]) + \
".mates." + outfmt, "wb") # mate output fasta file
        singouthandle = open(str(os.path.splitext(outfile)[0]) + \
".singletons." + outfmt, "wb") # singleton output fasta file
maphandle = open(mapfile, "rU") # mapping file, tsv format
if not suppress:
    uahandle = open(str(os.path.splitext(outfile)[0]) + \
    ".unassigned." + outfmt, "wb") # unassigned output fasta file



if verbose :
    print "\nDemultiplexing run started " + strftime("%Y-%m-%d %H:%M:%S") + "."
    if paired:
        print "Building dictionary of reads to process " + \
strftime("%Y-%m-%d %H:%M:%S") + "..."

# Set up iterators, and if paired, identify reads and positions where matches
# occur between the two files.
firstiter = SeqIO.parse(infile, infmt)
if paired:
    # Create two radix tries containing readnames from each file.
    #trie1=patricia()
    #for i,seq in enumerate(firstiter):
    #    trie1.addWord(getname(seq))
    #    if verbose:
    #        sys.stdout.write("Added from file 1: " + str(i+1))
    #        sys.stdout.flush()
    #        sys.stdout.write("\b"*len("Added from file 1: " + str(i+1)))
    #nmax=i+1
    #mateiter = SeqIO.parse(pairfile, infmt)
    #trie2=patricia()
    #print ""
    #for i,seq in enumerate(mateiter):
    #    seqname=getname(seq)
    #    if not trie1.isWord(seqname):
    #        nmax+=1
    #    trie2.addWord(seqname)
    #    if verbose:
    #        sys.stdout.write("Added from file 2: " + str(i+1))
    #        sys.stdout.flush()
    #        sys.stdout.write("\b"*len("Added from file 2: " + str(i+1)))
    #print ""
    
    # Create dictionary of names and positions.
    if verbose:
        print "Processing file 1..."
    readnames = dict((getname(seq, split_on, rn_part),[i,None]) for \
    i,seq in enumerate(firstiter))
    mateiter = SeqIO.parse(pairfile, infmt)
    if verbose:
        print "Processing file 2..."
    for i,seq in enumerate(mateiter):
        try:
            readnames[getname(seq, split_on, rn_part)][1]=i      
        except KeyError:
            readnames[getname(seq, split_on, rn_part)]=[None, i]
    nmax=len(readnames)
    iter_range = range(0,nmax)
    firstiter = SeqIO.parse(infile, infmt)#reload as exhausted by list comp.
    mateiter = SeqIO.parse(pairfile, infmt)#reload as exhausted by list comp.
    if verbose:
        print "Built " + strftime("%Y-%m-%d %H:%M:%S") + "."
        print str(nmax) + " mated and singleton reads to process."
else:
    if verbose:
        print "Counting reads in input file..."
    #get number of reads to process
    iter_range=firstiter
    #reinitialize generator (it will have been consumed in previous step)
    firstiter=SeqIO.parse(infile, infmt)
        
# Extract information from the mapping file.
if verbose:
    print "Parsing mapping file " + \
strftime("%Y-%m-%d %H:%M:%S") + "..." 
bcs = {} #Initialize bcs dictionary. Will have barcode as keyword 
#and position in list as value
ids = [] #Initialize list of sample IDs
primers = [] #Initialize list of primers
revprimers = []
line_ct = -1
#Extract barcodes, primers and IDs from mapping file
for line in maphandle :
    if line[0] != "#" :
        line_ct += 1
        separated_raw = line.split('\t')
        separated = [piece.rstrip() for piece in separated_raw]
        ids.append(separated[0])
        bcs[separated[1]] = line_ct
        try:
            primers.append(separated[2])
        except IndexError:
            pass
        if paired:
            #test for presence of data in the expected primer column
            try:
                #test that data is an ambiguous DNA sequence
                for nuc in separated[3]:
                    if nuc in ["A","C","G","T", "K", "M", "R", "Y", "S", "W",
                               "B", "V", "H", "D", "X", "N"]:
                        pass
                    else:
                        #if not DNA, fill with None
                        revprimers.append(None)
                        break
                revprimers.append(separated[3])
            except IndexError:
                pass
#Test that mapping file yielded some data.

try :
    ids[0]
except IndexError:
    print "No data loaded from mapping file. Check and try again."
    sys.exit(1)
if verbose:
    print "Parsed " + strftime("%Y-%m-%d %H:%M:%S") + "." 
#construct radix trie containing (mismatched) barcodes for efficient
#storage and searching
if verbose:
    print "Building radix trie from DNA barcodes " + strftime("%Y-%m-%d %H:%M:%S") + "..."
bc_trie = mmDNAtrie(bcs, max_bc_mismatch) # from demultiplex
if verbose:
    print "Built " + strftime("%Y-%m-%d %H:%M:%S") + "."
#Initializing counters:

# id_counts = matrix storing number of reads in each ID from mapping file [0]
id_counts = zeros( (len(ids), 2), dtype=object)
for i, id in enumerate(ids) :
    id_counts[i, 0] = id


# count_u = unassigned reads
count_u = 0
# count_a = assigned reads [initial value for output formatting of number of reads in category]
count_a = args.start_numbering_at[0]-1
# count_s = count of singleton reads, i.e., those that were assigned but had no mates
count_s = 0
# count_pm = sequences with primer mismatches,
#including unassigned that were rejected for too many primer mismatches
count_pm = 0
#count_norev = assigned reverse reads where no reverse primer was found when
#one had been specified in the maping file. May indicate a chimeric sequence.
count_norev = 0
i=0

if verbose:
    print "Demultiplexing..."
# Beginning of run, if so, get first two position values.
try:
    seq1 = firstiter.next()
except StopIteration:
    print "Error in input (-i) file. Please check and run again."
    sys.exit(2)
if paired:
    try:
        seq2 = mateiter.next()
    except StopIteration:
        print "Error in mate (-m) file. Please check and try again."
        sys.exit(3)
# Run identification loop.
for s in iter_range:
    i+=1
    if paired:
        # Determine which records to demultiplex.
        name1=getname(seq1, split_on, rn_part)
        name2=getname(seq2, split_on, rn_part)
        if name1==name2: # If current names are the same, then they are a pair.
            # Assign current seq1 to "record" to be identified.
            record = seq1
            # Assign current seq2 to "materecord" to be identified.
            materecord = seq2
            # Set next values to be tested in the following iteration.
            try:
                seq1 = firstiter.next()
                seq2 = mateiter.next()
            except StopIteration:
                pass
        else: # If the two current names are different, need to determine which
              # file to iterate through until one is found that is equal.
            if readnames[name1][1]==None:# If seq1 has no match in mate file,
                                         # iterate through it.
            #if not trie2.isWord(name1):
                record=seq1
                materecord=None
                try:
                    seq1=firstiter.next()
                except StopIteration:
                    pass
            elif readnames[name2][0]==None:# If seq2 has no match in first file
                                           # iterate through it.
            #if not trie1.isWord(name2):
                record=seq2
                materecord=None
                try:
                    seq2=mateiter.next()
                except StopIteration:
                    pass
    else:# If not paired, run through file.
        try:
            record = firstiter.next()
        except StopIteration or IndexError:
            break    
    #########################
    # Attempt to identify sequence stored in "record"
    result = identify_read(bc_trie, record, bcs, ids, primers, \
            max_pos=max_bc_st_pos, max_p_mismatch=max_prim_mismatch, \
            rpad=right_padding, bc_len=barcode_length, bc_type=barcode_type)
    if result[1] != "Unassigned" :#if a sample ID was assigned:
        if max_length:#set position of last nucleotide in read, based on
            #max_length parameter, if set.
            max_trim_pos = result[0]+max_length
        else:
            max_trim_pos = None
        #assign this record as forward read
        fwdseq = record[result[0]:max_trim_pos]
        fwdseq.description = "%s orig_bc=%s new_bc=%s bc_diffs=%i" \
                        % (fwdseq.name, result[2], result[3], \
                        ambiguous_seq_dist(result[2], result[3]))                    
        fwdseq.id = "%s_%i" % (result[1], count_a)
        if not paired:
            #write to outhandle if a
            #perfroming single-end run
            id_counts[bcs[result[3]],1] += 1
            count_a += 1
            SeqIO.write(fwdseq, outhandle, outfmt)
        elif materecord:
            #if paired-end run and materecord found,
            #write fwdseq to outhandle and scan revseq to trim primers
            revseq = materecord
            revseq.description = revseq.name
            revseq.id = "%s_%i" % (result[1], count_a)#id from matched fwdseq
            try:
                query = revprimers[bcs[result[3]]] #primer to look for.
                    #if a primer is present, it will search for it.
                try :
                    spos, epos = ambiguous_search(revseq.seq, query, max_prim_mismatch)
                    if max_length:
                        max_trim_pos = epos+max_length
                    else:
                        max_trim_pos = None
                    revseq = revseq[epos:max_trim_pos]
                #if no primer, it will leave the revseq unmodified and write it
                except TypeError:#returned if no match.
                    revseq.description = "%s %s" % (revseq.name, "rev. primer not found")
                    count_norev+=1                          
            except IndexError:#returned if no reverse primer.
                pass
            id_counts[bcs[result[3]],1] += 1
            count_a += 1
            SeqIO.write(fwdseq, outhandle, outfmt)
            SeqIO.write(revseq, mateouthandle, outfmt)
        else:
            #if paired-end run and no materecord, it will write to the
            #singletons file
            count_s += 1
            SeqIO.write(fwdseq, singouthandle, outfmt)
    #if no sample ID could be assigned to this read, scan paired mate
    #if using paired-ends, otherwise write to unassigned file.
    else:
        #if unassigned due to primer mismatches, increment relvant counter
        #result[4] will only be full when unassigned if there were too many
        #primer mismatches
        if result[4]:
            count_pm += 1
        #if not paired, means read is unassigned and can be written to
        #the file of unassigned reads
        if not paired:
            count_u += 1
            if not suppress:
                SeqIO.write(record, uahandle, outfmt)
        #if it is paired, but the mate record is not present,
        #the read can be written to the unassigned file as it's not
        #possible to identify where it came from.
        elif paired and not materecord:
            count_u += 1
            if not suppress:
                SeqIO.write(record, uahandle, outfmt)
        #if it is paired and a materecord exists, the materecord
        #can be scanned to try and identify the read's origin
        else:
            result = identify_read(bc_trie, materecord, bcs, ids, primers, \
            max_pos=max_bc_st_pos, max_p_mismatch=max_prim_mismatch, \
            rpad=right_padding, bc_len=barcode_length, bc_type=barcode_type)
            #if it gets assigned an id, write it to the outhandle and use
            #the new id (fwdseq.id) with the reverse read
            if result[1] != "Unassigned" :
                if max_length:#set position of last nucleotide in read,
                    #based on max_length parameter, if set.
                    max_trim_pos = result[0]+max_length
                else:
                    max_trim_pos = None
                #assign this record as forward read
                fwdseq = materecord[result[0]:max_trim_pos]
                fwdseq.description = "%s orig_bc=%s new_bc=%s bc_diffs=%i" \
                                % (fwdseq.name, result[2], result[3], \
                                ambiguous_seq_dist(result[2], result[3]))                    
                fwdseq.id = "%s_%i" % (result[1], count_a)
                revseq = record
                revseq.description = revseq.name
                revseq.id = "%s_%i" % (result[1], count_a)#id from matched fwdseq
                try:
                    query = revprimers[bcs[result[3]]] #primer to look for
                    #if a primer is present, it will search for it
                    try :
                        spos, epos = ambiguous_search(revseq.seq, query, max_prim_mismatch)
                        if max_length:
                            max_trim_pos = epos+max_length
                        else:
                            max_trim_pos = None
                        revseq = revseq[epos-1:max_trim_pos]
                    #if no primer, it will leave the revseq unmodified and write it
                    except TypeError:
                        revseq.description = "%s %s" % (revseq.name, "rev. primer not found")
                        count_norev+=1
                except IndexError:
                    pass
                id_counts[bcs[result[3]],1] += 1
                count_a += 1
                SeqIO.write(fwdseq, outhandle, outfmt)
                SeqIO.write(revseq, mateouthandle, outfmt)
            #if no match was found, can write both records to the unassigned file.
            else:
                count_u += 1
                if not suppress:
                    SeqIO.write(record, uahandle, outfmt)
                    SeqIO.write(materecord, uahandle, outfmt)
                   
#Progress indicator
    if verbose:
        sys.stdout.write("Assigned: " + str(count_a) + ", Processed: " + str(i+1))
        sys.stdout.flush()
        sys.stdout.write("\b"*len("Assigned: " + str(count_a) + ", Processed: " + str(i+1)))
#set total count after run finished
count = i+1



logpath = open(str(os.path.splitext(outfile)[0]) + ".log","wb")
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
"Unassigned reads: " + str(count_u) + "\n")
if paired:
    logpath.write("Singletons: " + str(count_s) + "\n"
                  "Reverse reads without reverse primer (may indicate chimeras): "
                  + str(count_norev) + "\n\n"
                  "Sample id.\tPaired Reads\n")
else:
    logpath.write("Sample id.\tReads\n")
savetxt(logpath, id_counts, fmt='%s', delimiter='\t')


if paired:
    mateouthandle.close()
    singouthandle.close()
        
outhandle.close()
maphandle.close()
if not suppress:
    uahandle.close()
logpath.close()
if verbose:
    print "\nRun finished " + strftime("%Y-%m-%d %H:%M:%S") + "."





