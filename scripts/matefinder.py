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
    description='''Takes a demultiplexed next-gen sequencing fasta file for one
    direction in a paired-end run and demultiplexes the other direction's
    sequence file by matching the original read names.''',
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
matedata = SeqIO.index(matefile, matefmt)

print "Finding mates..."
mated_ctr = 0
for i, rec in enumerate(indata):
    if verbose :#print line number to stdout
        sys.stdout.write("Processing read: " + str(i+1))
        sys.stdout.flush()
        sys.stdout.write("\b"*len("Processing read: " + str(i+1)))
    #get original read name and direction identifier which should be the second
    #part of the post-demultiplexing id.
    split_id = rec.description.split()
    orig_id, ending = split_id[1].split("/")
    #depending on the direction identifier ("ending") of the demultiplexed read,
    #construct the full original read name of the mate
    if int(ending) == dir_ids[0]:
        mate_orig_id = '/'.join([orig_id, str(dir_ids[1])])
    elif int(ending) == dir_ids[1]:
        mate_orig_id = '/'.join([orig_id, str(dir_ids[0])])
    else:
        print "Read direction identifier doesn't match any that were entered"
        sys.exit(1)
    #look for the original mate id in the matedata keys and then set it as new
    #sequence to be written, incorporating the appropriate modifications to its
    #length and header.
    try:
        new_rec = matedata[mate_orig_id]
        #if mapping file specified containing primers to trim:
        if mapfile:
            for query in trimlist:
                try :
                    spos, epos = ambiguous_search(new_rec.seq, query, max_mismatch)
                    if max_length:
                        max_trim_pos = epos+max_length
                    else:
                        max_trim_pos = None
                    new_rec = new_rec[epos:max_trim_pos]
                    break
                except TypeError:
                    pass        
        #if no mapping file, set the trimming based on the "trim"
        #and "max_length" parameters, if present.
        else:
            if max_length:
                if trim:    
                    max_trim_pos = trim+max_length
                else:
                    max_trim_pos = max_length
            else:
                max_trim_pos = None
            new_rec = new_rec[trim:max_trim_pos]
        #create new headers for mate sequences from demultiplexed names and then
        #write to output file
        new_rec.id = "%s" % split_id[0]
        new_rec.description = mate_orig_id
        SeqIO.write(new_rec, outhandle, outfmt)
        mated_ctr += 1
    except KeyError:
        pass
total_reads = i+1
orphans = i+1-mated_ctr


# Write log
logpath = open(str(os.path.splitext(outfile)[0]) + ".log","wb")
logpath.write("Logfile for mate-finding of " \
+ infile + "and" + matefile + "\n" + strftime("%Y-%m-%d %H:%M:%S") + "\n\n" \
"Parameters specified:\n" \
"Mapping file: " + str(mapfile) + "\n" \
"Maximum mismatches with sequences to be trimmed: " + str(max_mismatch) + "\n" \
"Nucleotides to trim from left: " + str(trim) + "\n" \
"Maximum length: " + str(max_length) + "\n\n" \
"Counts:"
"\nTotal reads parsed: " + str(total_reads) + "\n" \
"Reads mated: " + str(mated_ctr) + "\n" \
"Orphan reads:" + str(orphans) + "\n")
logpath.close()

print "\n\nLog file written (" + str(os.path.splitext(outfile)[0]) + ".log" + ")\n"
    
#inhandle.close()
#matehandle.close()
try:
    maphandle.close()
except NameError:
    pass    
outhandle.close()
logpath.close()
    
print "Run finished " + strftime("%Y-%m-%d %H:%M:%S") + "."

