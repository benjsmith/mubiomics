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
from collections import deque, Counter

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
#parser.add_argument('-d', '--direction_ids', required=False, nargs='+',
#                    default='',
#                    type=int, help='''enter numbers in read names identifying
#                    direction of read. For Illumina HiSeq, this is the last
#                    number after the "/" in the read names e.g.,usually single
#                    end runs use "/0", so you would specify "-d 0".
#                    If demultiplexing a paired-end run, enter two numbers
#                    separated by a space, e.g., if forward reads end in "/1"
#                    and reverse reads end in "/3", you would enter "-d 1 3".
#                    For Illumina MiSeq, the read direction is the number
#                    before the ":N" towards the end of the read name. For 454,
#                    there is no need to set this parameter.''')
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
#parser.add_argument('-u', '--use_indexdb', action='store_true',
#                    help='''activate database indexing mode. This should be
#                    used if the combined size of the input files on a paired run
#                    is larger than two thirds the available RAM. It parses the
#                    entries and stores them in an SQL database file in the directory
#                    where the input sequences reside. This reduces memory
#                    consumption, but also significantly reduces the speed of the
#                    run, especially on slow hard drives. For single end runs, this
#                   is never necessary.''')
#parser.add_argument('-x', '--index_exists', action='store_true',
#                    help='''activate this if sequences have already been indexed
#                    and stored in an SQL database file by this script using
#                    index_db, or by another script using Biopython's
#                    SeqIO.index_db(...). If this is activated, make your input
#                    file (-i/--input_fp) the path to the database file, there
#                    is no need to activate -u/--use_indexdb in this case, but
#                    you should still specify the format of the sequences in the
#                    index file using the -I/--in_fmt option.''')
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
        print "--barcode_type (-b) specified is an unsupported type"
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
        print "Building list of readnames to process " + \
strftime("%Y-%m-%d %H:%M:%S") + "..."

def idx(myitem, mylist):
    '''Index function that returns None when item is not present.
    
    '''
    try:
        pos = mylist.index(myitem)
        mylist.pop(pos)
        return pos
    except ValueError:
        return None
    
def getname(seqrecord, splitchar=split_on, partn=rn_part):
    '''Extracts the name of the read based on sequencing platform.
    
    '''
    return seqrecord.description.split(splitchar)[partn]
    

## New Approach ##
firstiter = SeqIO.parse(infile, infmt)
firstreadnames = [getname(seq, split_on, rn_part) for seq in firstiter]
if not paired:
    nmax=len(firstreadnames)
    if verbose:
        print "Built " + strftime("%Y-%m-%d %H:%M:%S") + "."
        print str(nmax) + " reads to process."
else:
    mateiter = SeqIO.parse(pairfile, infmt)
    matereadnames = [getname(seq, split_on, rn_part) for seq in mateiter]
    allreadnames=firstreadnames[:]
    allreadnames.extend(matereadnames)
    nmax=len(set(allreadnames))
    allreadnames=[]
    firstiter = SeqIO.parse(infile, infmt)#reload as exhausted by list comp.
    mateiter = SeqIO.parse(pairfile, infmt)#reload as exhausted by list comp.
    if verbose:
        print "Built " + strftime("%Y-%m-%d %H:%M:%S") + "."
        print str(nmax) + " reads to process."
        #print "Getting positions " + strftime("%Y-%m-%d %H:%M:%S") + "..."
    ##position generators
    #posmateinfirst = (idx(name, firstreadnames) for name in matereadnames)
    #posfirstinmate = (idx(name, matereadnames) for name in firstreadnames)
    # position lists generated by map
    #posmateinfirst = map((lambda name: idx(name, firstreadnames)), matereadnames)
    #if verbose:
    #    print "First set of positions obtained " + \
#strftime("%Y-%m-%d %H:%M:%S") + "." 
    #posfirstinmate = map((lambda name: idx(name, matereadnames)), firstreadnames)
    #clear readname lists if using position lists
    #firstreadnames=[]
    #matereadnames=[]
    #if verbose:
    #    print "All positions obtained " + \
#strftime("%Y-%m-%d %H:%M:%S") + "." 
    

    
        
#if paired:
#    #make set of sequence base names that are in both files by stripping directional
#    #identifier
#    readnames  = []
#    rnappend=readnames.append
#    dir_ids = []
#    for i, n in enumerate(SeqIO.parse(infile, infmt)):
#        if verbose:
#            sys.stdout.write("File 1, Read: " + str(i+1))
#            sys.stdout.flush()
#            sys.stdout.write("\b"*len("File 1, Read: " + str(i+1)))
#        seqname1 = n.description.split(split_on)[rn_part]
#        #readnames1.append(seqname1)
#        rnappend(seqname1)
#        # In case there are more characters after the read number e.g. for MiSeq,
#        # extract and store them to build up the list of read direction ids.
#        ending = split_on.join(n.description.split(split_on)[rn_part+1:])
#        if ending not in dir_ids:
#            dir_ids.append(ending)
#    print ""
#    for j, m in enumerate(SeqIO.parse(pairfile, infmt)):
#        if verbose:
#            sys.stdout.write("File 2, Read: " + str(j+1))
#            sys.stdout.flush()
#            sys.stdout.write("\b"*len("File 2, Read: " + str(j+1)))
#        seqname2 = m.description.split(split_on)[rn_part]
#        rnappend(seqname2)
#        ending = split_on.join(m.description.split(split_on)[rn_part+1:])
#        if ending not in dir_ids:
#            dir_ids.append(ending)
#    # For index(_db), we need to define a function to create paired sets of
#    # dict keys since these methods only take the id, which needs to be
#    # supplemented with read direction information. For e.g. MiSeq
#    # sequence files, the id doesn't contain any info about the read
#    # direction, however, we can get this from the dir_ids we stored earlier.
#    readCtr=Counter(readnames) #used in makeid function
#    readnames=readCtr.keys() #used in makeid function
#    midx=readnames.index # define function so it's not revaluated on each
#    # makeid call during sequence indexing.
#    def makeid(name, idx_fn=midx, dir_ids=dir_ids): # make local variables for speed.
#        pos=idx_fn(name) # idx_fn defined after readnames is built
#        cnt=readCtr[name]
#        val=logiclist[pos]
#        if cnt==2:
#            if val==0:
#                logiclist[pos]=1
#            else:
#                logiclist[pos]=0
#        id = split_on.join([name, dir_ids[val]])
#        return id
#    
#    if verbose:
#        print ""
#        print "Built " + strftime("%Y-%m-%d %H:%M:%S") + "."
#    # Make logiclist form makeid function
#    logiclist=[0]*len(readnames)    
#    if not idx_exists:
#    #if paired, index sequences
#        if verbose:
#            print "Indexing input sequence files..."
#        # index sequence input files using SQL database, creating
#        # new index file.
#        indexfile = str(os.path.splitext(os.path.abspath(infile))[0]) + ".idx"
#        try:
#            os.remove(indexfile)
#        except OSError:
#            pass
#    # index sequences using SQL db or in memory
#    if use_indexdb:
#        indata = SeqIO.index_db(indexfile, [infile, pairfile], infmt,
#                                key_function=makeid)
#    else:
#        initindata = SeqIO.index(infile, infmt,
#                                 key_function=makeid)
#        pairdata = SeqIO.index(pairfile, infmt,
#                               key_function=makeid)
#        indata = MultiIndexDict(initindata, pairdata)
##elif idx_exists and in_mem:
#    #    indata = SeqIO.index_db(":memory:", infile)
#    #else: #if index SQL db already exists, read it in.
#    #    # Make logiclist for makeid function.
#    #    logiclist=[0]*len(readnames)
#    #    indata = SeqIO.index_db(infile, key_function=makeid)
## if not paired, will parse the file on the fly further down the script.

if verbose:
    print "Parsing mapping file " + \
strftime("%Y-%m-%d %H:%M:%S") + "..." 
#logiclist=[0]*len(readnames) #Reset logiclist so seq index __getitem_ calls
# later on are correctly synchronized to the number of ends each read has.
# The __getitem__ method calls the key_function, makeid
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
        separated = line.split()
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
    print "Building radix trie " + strftime("%Y-%m-%d %H:%M:%S") + "..."
bc_trie = mmDNAtrie(bcs, max_bc_mismatch) # from demultiplex
if verbose:
    print "Built " + strftime("%Y-%m-%d %H:%M:%S") + "."

#Initializing counters:
# id_counts = matrix storing number of reads in each ID from mapping file [0]
# count_u = unassigned reads
# count_a = assigned reads [initial value for output formatting of number of reads in category]
# count_s = count of singleton reads, i.e., those that were assigned but had no mates
# count_pm = sequences with primer mismatches, including unassigned that were rejected for too many primer mismatches

id_counts = zeros( (len(ids), 2), dtype=object)
for i, id in enumerate(ids) :
    id_counts[i, 0] = id

count_bad = 0
count_u = 0
count_a = args.start_numbering_at[0]-1
count_s = 0
count_pm = 0

if verbose:
    print "Demultiplexing..."
# 
r = 0
#if not paired:
#    readnames = SeqIO.parse(infile, infmt)
#while readnames:
#    if paired:
#        try: #try to take first element of readnames
#            recname = readnames[readname_pos]
#        except IndexError: #if it is exhausted, break.
#            break
#        try :
#            #if paired-end, the selected record may have the other dir_id
#            record = indata[split_on.join([recname, str(dir_ids[0])])]
#            try:
#                #print split_on.join([recname, str(dir_ids[1])])
#                materecord = indata[split_on.join([recname, str(dir_ids[1])])]
#            except KeyError:
#                materecord = None
#                #mate pair is absent
#        except KeyError:
#            try:
#                #if get here, the mate pair is absent
#                record = indata[split_on.join([recname, str(dir_ids[1])])]
#                materecord = None
#            except KeyError:
#                #if get here, wrong read direction identifiers were entered.
#                print "Read direction identifiers in reads (-d/--direction_ids) \
#don't match any that were entered. Run 'demultiplexer.py -h' for help."
#                sys.exit(1)
#        readname_pos+=1
#    else:
#        try:
#            record = readnames.next()
#        except StopIteration:
#            break
#    # A record is always assigned successfully
    
## New Approach ##
pos1_idx=0
pos2_idx=0
for i in range(0,nmax):
    if paired:
        # Get positions.
        # Test if beginning of run, if so, get first two position values.
        if pos1_idx==0 and pos2_idx==0:
            seq1 = firstiter.next()
            seq2 = mateiter.next()
            pos1 = idx(getname(seq1), matereadnames)
            pos1_idx+=1
            pos2 = idx(getname(seq2), firstreadnames)
            pos2_idx+=1
        # Determine which records to demultiplex.
        # If both current positions match both current indexes, assign both
        # reads.
        if pos1==pos2_idx and pos2==pos1_idx:
            try:
                record = seq1
                seq1 = firstiter.next()
                pos1 = idx(getname(seq1), matereadnames)
                pos1_idx+=1
            except StopIteration or IndexError:
                pass
            try:
                materecord = seq2
                seq2 = mateiter.next()
                pos2 = idx(getname(seq2), firstreadnames)
                pos2_idx+=1
            except StopIteration or IndexError:
                pass
        elif pos1_idx<pos2 or pos1==None:
            # Maintain previous values of seq2 and pos2.
            # Increment seq1 and pos1.
            try:
                record = seq1
                seq1 = firstiter.next()
                pos1 = idx(getname(seq1), matereadnames)
                pos1_idx+=1
                pos2_idx-=1
                materecord = None
            except StopIteration or IndexError:
                pass
        elif pos1>pos2_idx or pos2==None:
            # Maintain previous values of seq1 and pos1.
            # Increment seq2 and pos2.
            try:
                #increment position for next iteration
                record = seq2
                seq2 = mateiter.next()
                pos2 = idx(getname(seq2), firstreadnames)
                pos2_idx+=1
                pos1_idx-=1
                materecord = None
            except StopIteration or IndexError:
                pass
        else:
            print "Paired run exited because only one of the position \
indices matched the current position, where both indices should match both \
or one current position should be None or greater than the other position index."
            print "pos1: " + str(pos1)
            print "pos1_idx: " + str(pos1_idx)
            print "pos2: " + str(pos2)
            print "pos2_idx: " + str(pos2_idx)
            sys.exit(1)
        
    else:
        try:
            record = firstiter.next()
        except StopIteration or IndexError:
            break    
    #########################
    # Attempt to identify sequence stored in "record"
    result = identify_read(bc_trie, record, bcs, ids, primers, \
            max_pos=max_bc_st_pos, max_p_mismatch=max_prim_mismatch, \
            rpad=right_padding, bc_len=barcode_length, bc_type=barcode_type)
    #print materecord
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
            query = revprimers[bcs[result[3]]] #primer to look for
            #if a primer is present, it will search for it
            try :
                spos, epos = ambiguous_search(revseq.seq, query, max_prim_mismatch)
                if max_length:
                    max_trim_pos = epos+max_length
                else:
                    max_trim_pos = None
                revseq = revseq[epos:max_trim_pos]
            #if no primer, it will leave the revseq unmodified and write it
            except TypeError:
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
                   
    r += 1
#Progress indicator
    if verbose:
        sys.stdout.write("Assigned: " + str(count_a) + ", Processed: " + str(r))
        sys.stdout.flush()
        sys.stdout.write("\b"*len("Assigned: " + str(count_a) + ", Processed: " + str(r)))
#set total count after run finished
count = r



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
    logpath.write("Singletons: " + str(count_s) + "\n\n"
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





