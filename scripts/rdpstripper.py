#!/usr/bin/env python

#Author: Dr. Benjamin C. Smith, 2011

import sys, os, re, argparse, csv, fileinput
from collections import deque, Counter
import cPickle as Pickle
from numpy import zeros, hstack, vstack, array, arange




parser = argparse.ArgumentParser(
    description='''Convert a csv file produced by RDP multiclassifier, with
"fixRank" format, to an OTU table .''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-i','--input', nargs=1, required=True,
                    type=str, help='''input file. Should be a csv file
                    containing a list of classifications produced by
                    an sqlite3 command on a guppy classify database. To
                    get the sqlite3 commands, run script with parameter
                    -c/--command.''')
parser.add_argument('-o', '--outpath', nargs='?', type=str,
                    help='''output path. The resulting OTU table is written
                    to this filename. If not set, creates a file in the current
                    directory called results_table.txt''')
parser.add_argument('-r', '--rank', nargs=1, type=str, default=['genus'],
                    help='''specify taxonomic rank at which to extract counts.
                    ''')
parser.add_argument('-c', '--conf', nargs=1, default=[0.8], type=float,
                    help='''set minimum acceptable bootstrap confidence for
                    an assignment to be counted.''')
parser.add_argument('-s', '--split_on', nargs=1, default=['_'], type=str,
                    help='''set the chracter to split the sample name on. E.g.,
                    if the read name is sample1_00001 and the split is set to
                    "_", the sample name will be read as "sample1".''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='''activate verbose mode. Indicates progress and other
                    information.''')

#Setting required variables from command line arguments.
args = parser.parse_args()
infile = args.input[0]
outpath = args.outpath
rank = args.rank[0]
conf = args.conf[0]
split_on = args.split_on[0]
verbose = args.verbose


if not outpath:
    outpath = '/'.join(os.path.abspath(infile).split('/')[:-1]) + \
    "/results_table.txt"
try:
    os.mkdir('/'.join(outpath.split('/')[:-1]))
except OSError:
    pass

print "Running..."
    
inhandle = csv.reader(open(infile, "rU"), delimiter='\t')#input file handle
outhandle = csv.writer(open(outpath, "wb"), dialect='excel', delimiter='\t')# output file handle

#generate dictionary with keys = sample names, values = list of otus in sample
samples = {}
for j, row in enumerate(inhandle) :
    if verbose :#print line number to stdout
        sys.stdout.write("Read: " + str(j))
        sys.stdout.flush()
        sys.stdout.write("\b"*len("Read: " + str(j)))
    sample = str(row[0].split(split_on)[0])
    for i, part in enumerate(row):#select only columns pertaining to specified rank
        if part==rank:
            col=i
            if row[col+1]>=conf:#test that the confidence is high enough
                tax_name = row[col-1]
    try:
        samples[sample].append(tax_name)
    except KeyError:
        samples[sample] = [tax_name]

#generate list of sample ids
#generate list of all observed otus (no duplicates)
sample_ids = []
otus = []
for key, otu_list in samples.items():
    sample_ids.append(str(key))
    for otu in otu_list:
        if otu not in otus:
            otus.append(otu)

#create empty count matrix
count_matrix = zeros( (len(otus), len(sample_ids)), dtype=int )

#populate count matrix with counts from samples dictionary
for key, values in samples.items():
    col = sample_ids.index(key)
    for value in values:
        row = otus.index(value)
        count_matrix[row, col] += 1

#construct output table
otu_ids = ["#OTU ID"] + list(arange(len(otus)))
otus = ["Consensus Lineage"] + otus
sample_ids = array(sample_ids)
len_otus = len(otus)
otus = array([[i] for i in otus])
otu_ids = array([[i] for i in otu_ids])
pre_table = vstack((sample_ids, count_matrix))
final_table = hstack((otu_ids, pre_table, otus))

outhandle.writerow(["#rdp multiclassifier OTU table"])
for row in final_table:
    outhandle.writerow(row)
            

print "\nFinished."

