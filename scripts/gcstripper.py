#!/usr/bin/env python

#Author: Dr. Benjamin C. Smith, 2011

import sys, os, re, argparse, csv, fileinput
from collections import deque, Counter
import cPickle as Pickle
from numpy import zeros, hstack, vstack, array, arange




parser = argparse.ArgumentParser(
    description='''Convert a csv file produced by SQL queries on a guppy classify
    sqlite3 databse, to an OTU table .''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-i','--input', nargs='?',
                    type=str, help='''input file. Should be a csv file
                    containing a list of classifications produced by
                    an sqlite3 command on a guppy classify database. To
                    get the sqlite3 commands, run script with parameter
                    -c/--command.''')
parser.add_argument('-o', '--outpath', nargs='?', type=str,
                    help='''output path. The resulting OTU table is written
                    to this filename. If not set, creates a file in the current
                    directory called results_table.txt''')
parser.add_argument('-c', '--command', action='store_true',
                    help='''display sqlite3 command that was used to
                    generate .''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='''activate verbose mode. Indicates progress and other
                    information.''')

#Setting required variables from command line arguments.
args = parser.parse_args()
infile = args.input
outpath = args.outpath
verbose = args.verbose
command = args.command

if command:
    print"\n#################"
    print "sqlite 3 command for generating a csv file containing \
classifications."
    print "#################"
    print "The variables\nDB=<path_to_database>\n\
THRESHOLD=<minimum_required_likelihood> \n\
TAX=<required_taxonomic_level>\n\
must be set prior to running this command."
    print "################\n"
    print 'sqlite3 -header -csv "${DB}" " \n\
SELECT pc.placement_id, \n\
       placement_names.origin, \n\
       taxa.tax_name, \n\
       pc.desired_rank, \n\
       pc.likelihood \n\
FROM placement_classifications AS pc \n\
INNER JOIN taxa \n\
ON pc.tax_id=taxa.tax_id \n\
INNER JOIN placement_names \n\
ON pc.placement_id=placement_names.placement_id \n\
WHERE pc.likelihood>"${THRESHOLD}" \n\
AND pc.desired_rank='"${TAX}"' \n\
AND pc.rank='"${TAX}"' \n\
ORDER BY pc.placement_id \n\
" >"${TAX}"_at_"${THRESHOLD}".csv\n'
    sys.exit(0)
elif not input:
    print "usage: gc_stripper.py [-h] -i INPUT/-c [-o [OUTPATH]] [-v]"
    print "gc_stripper.py: error: if argument -c not given, \
argument -i/--input is required"
    sys.exit(1)

if not outpath:
    outpath = '/'.join(os.path.abspath(infile).split('/')[:-1]) + \
    "/results_table.txt"
try:
    os.mkdir('/'.join(outpath.split('/')[:-1]))
except OSError:
    pass

print "Running..."
    
inhandle = csv.reader(open(infile, "rU"), delimiter=',')#input file handle
outhandle = csv.writer(open(outpath, "wb"), dialect='excel', delimiter='\t')# output file handle



inhandle.next()#skip first row of input file containing headers
#generate dictionary with keys = sample names, values = list of otus in sample
samples = {}
for j, row in enumerate(inhandle) :
    if verbose :#print line number to stdout
        sys.stdout.write("Read: " + str(j))
        sys.stdout.flush()
        sys.stdout.write("\b"*len("Read: " + str(j)))
    sample = str(row[1])
    tax_name = str(row[2])
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

outhandle.writerow(["#guppy classify OTU table"])
for row in final_table:
    outhandle.writerow(row)
            

print "\nFinished."

