#!/usr/bin/env python
#Pools assigned OTUs with identical names and renumbers the remaining distinct
#OTUs. Also allows filtering out OTUs with less than "min_cts" in at least
#one sample.
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


import sys, os, re, argparse, csv
from numpy import array, sum, append, amax, hstack, savetxt, linspace
from itertools import product
from time import strftime


parser = argparse.ArgumentParser(description='''Filter an OTU table by pooling
OTUs with identical names. (Optionally: discard OTUs with less than a specified
minimum proportion of counts in any one sample)''')
parser.add_argument('-i','--infile', required=True, nargs=1,
		    type=str, help='''input filepath. Should be a tab-delimited
		    OTU table.''')
parser.add_argument('-o', '--outfile', required=True, nargs=1, type=str,
		    help='''output filepath. The resulting pooled OTU table is
		    written here.''')
parser.add_argument('-m', '--min', nargs=1, default=[0], type=float,
		    help='''set the lowest acceptable read proportion e.g. 0.01
		    will keep all reads that comprise more than 1 percent of a
		    sample's composition and set any below 1 percent to 0. Any
		    rows (OTUs) that conatin all 0s after filtering will be
		    discarded. ''')
parser.add_argument('-r', '--reads', action='store_true', 
                    help='''print information about number of reads''')
args = parser.parse_args()
min_cts = args.min[0]
if min_cts > 1 or min_cts < 0:
    print "Invalid minimum count threshold (-m/--min) set. \
Must be in range [0, 1]."
    sys.exit(1)
infile = args.infile[0]
outfile = args.outfile[0]


print "\nRun started " + strftime("%Y-%m-%d %H:%M:%S") + "."

inhandle = csv.reader(open(infile, 'rU'), delimiter='\t')
outhandle = csv.writer(open(outfile, 'wb'), delimiter='\t')
#Write header lines of output file.
outhandle.writerow([''.join(inhandle.next()) + ', pooled'])
outhandle.writerow(inhandle.next())

#collect sample names, using first line of file
inhandle = csv.reader(open(infile, 'rU'), delimiter='\t')
inhandle.next()
sample_ids = [column for column in inhandle.next() if \
re.search(column, "#OTU ID"'|'"Consensus Lineage")==None]
otu_names = []
otu_dict = {}

#build list of OTU names
inhandle = csv.reader(open(infile, 'rU'), delimiter='\t')
inhandle.next()# skip first line
inhandle.next()#skip second line
for line in inhandle :
    if line[-1] not in otu_names :
	otu_names.append(line[-1])
	# K,V = name of taxon, list of number of occurrences in each sample
	#there may be more than one V for each K.
	otu_dict[line[-1]] = [line[1:-1]]
    else :
	otu_dict[line[-1]].append(line[1:-1])

#create array of total counts per sample by summing along columns
total_counts_per_otu=array([item for sublist in otu_dict.values() \
			    for item in sublist], dtype=int).sum(axis=0)
min_acceptable=min_cts*total_counts_per_otu


otu_table = [] #empty array that will be filled with filtered count data
ictr = 1 #counter for assigning new OTU IDs.
#counters for tracking numbers of reads in intial and final OTU tables
tot_start_cts = 0
tot_end_cts = 0
for key in otu_dict.keys() :
#calculate per-sample counts for OTU name
    otu_counts = array(otu_dict[key], dtype=int).sum(axis=0)
    #total number of starting reads per row, before filtering
    tot_start_cts += otu_counts.sum()
    #replace with zero, any counts falling below the minimum acceptable
    #proportion
    otu_counts = array([[otu_val if otu_val>min_val else 0][0] for otu_val, \
	min_val in zip(otu_counts, min_acceptable)], dtype=int)
    #test that otu_counts is not a row of zeroes. If not, proceed
    if otu_counts.sum() > 0 :
        #create row for output 
        if key != 'Noise' : #if not the "Noise" OTU add otu_id from ictr
            # and increment it by 1.
            otu_id = array( [ictr], dtype=object)
            ictr += 1
        else: # if "Noise" OTU, set otu_id to '0' and don't increment ictr.
            otu_id = array( [0], dtype=object)
        otu_name = array( [key], dtype=object)
        otu_row = hstack( (otu_id, otu_counts, otu_name) )
        tot_end_cts += otu_counts.sum()
        otu_table.append(otu_row)
otu_table = array(otu_table) # convert to numpy array to allow easy sorting
final_otu_table = otu_table[otu_table[:,0].argsort(),:].tolist() # sort
#otu_table by otu_id and convert back to list
for row in final_otu_table:
    outhandle.writerow(row)
print "Finished.\n"
print "Final OTU table preview: "
print final_otu_table
 
# Write log
logpath = open(str(os.path.splitext(outfile)[0]) + ".log","wb")
logpath.write("Logfile for OTU pooling of " \
+ infile + "\n" + strftime("%Y-%m-%d %H:%M:%S") + "\n\n" \
"Parameters specified:\n" \
"Minimum read threshold: " + str(min_cts) + "\n" \
"Counts:"
"\nTotal reads in input OTU table: " + str(tot_start_cts) + "\n" \
"Total reads in output OTU table: " + str(tot_end_cts) + "\n" \
"Reads discarded by setting minimum threshold to " + str(100* min_cts) \
    + " %: " + str(tot_start_cts-tot_end_cts) + "\n")
logpath.close()

print "\n\nLog file written (" + str(os.path.splitext(outfile)[0]) + ".log" + ")\n"

if args.reads:
    print '\nTotal reads in input OTU table: ' + str(tot_start_cts)
    print 'Total reads in output OTU table: ' + str(tot_end_cts)
    print 'Reads discarded by setting minimum threshold to ' + str(100* min_cts) \
    + ' %: ' + str(tot_start_cts-tot_end_cts) + '\n'

