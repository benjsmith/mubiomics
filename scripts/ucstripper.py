#!/usr/bin/env python
#ucstripper.py
#convert results from a usearch database search to an OTU table for single-end
#data

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

import sys, os, re, argparse
from collections import defaultdict, deque
from itertools import product
from numpy import *


parser = argparse.ArgumentParser(
    description='''Strip data from a usearch .uc output file and produce an OTU
    table''',
    formatter_class= argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@')
parser.add_argument('-i','--input_fp', required=True, nargs=1,
                    type=str, help='''input filepath. Should be a .uc file
                    produced by a usearch --cluster run against a database.''')
parser.add_argument('-o', '--output_fp', required=True, nargs=1, type=str,
                    help='''output filepath. The resulting OTU table is written
                    here.''')
parser.add_argument('-m', '--map_fp', required=True, nargs=1, type=str,
                    help='''mapping filepath. Should be a tab-separated mapping
                    file with a sample name in the first column of each line.
                    The sample names in this file will be used as the column
                    headings in the output.''')
parser.add_argument('-s', '--split_on', nargs=1, default=['_'], type=str,
                    help='''set the chracter to split the sample name on. E.g.,
                    if the read name is sample1_00001 and the split is set to
                    "_", the sample name will be read as "sample1".''')
parser.add_argument('-M', '--min_id', nargs=1, default=[97.0], type=float,
                    help='''set minimum acceptable percent identity for a hit to
                    to be accepted.''')
parser.add_argument('-S', '--db_name_split', required=True, nargs=1,
                    default=["|"],
                    help='''set the separator character for reference sequence names.
                    This can be used for selecting only a part of the name of the
                    matched reference sequences as the taxon name.
                    E.g. If the name of a sequence is
                    "44919 Neisseria meningitidis str. M1976 AF310437.1 1..1471"
                    and you just wanted "Neisseria meningitidis", you could
                    set -S " " -P 2 3 (see notes for -P, below). If setting
                    -P "All" to take the whole name then you can ignore this
                    parameter.''')
parser.add_argument('-P', '--db_name_part', required=True, nargs='+',
                    default=["All"], type=int,
                    help='''set the part or parts to select from the
                    reference sequence name (see notes for -S, above). To
                    select all parts set -P "All". To select a range of parts
                    list the part numbers, separated by spaces. E.g., -P 1 2 3 to
                    select the 1st, 2nd and 3rd parts (numbering starts at 1).''')
parser.add_argument('-n', '--noise', action='store_true',
                    help='''condense counts for all clusters that do not match
                    any database sequences into a group called "Noise".''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='''activate verbose mode. Indicates progress and other
                    information.''')

#Setting required variables from command line arguments.
args = parser.parse_args()
infile = args.input_fp[0]
outfile = args.output_fp[0]
# tab-delimited meta-data file
mapfile = args.map_fp[0]
split_on = args.split_on[0]
min_id = args.min_id[0]
db_name_split=args.db_name_split[0]
db_name_part=args.db_name_part
if db_name_part!="All":
    db_name_part=[i-1 for i in db_name_part]
noise = args.noise
verbose = args.verbose


#open input/output handles
handle1 = open(infile, "rU") # uclust output file handle
handle2 = open(outfile, "wb") # output text file handle
handle3 = open(mapfile, 'rU') #mapping file handle

if verbose:
    print "Running..."

handle2.write("#ucstripper OTU Table\n") # add header line to outout file

#parse mapping file to build list of sample names, these will be the column
#headers for the count data.
sample_ids = ['#OTU ID']
for line in handle3 :
    if line[0] != "#" :
        id = line.split()[0].split(split_on)[0]
        if id not in sample_ids:
            sample_ids.append(id)
sample_ids.append('Consensus Lineage')

if verbose:
    print "Processing reads for " + str(len(sample_ids)-2) + " samples."

#parse the input, 'hit' by 'hit', and add info from columns to
#data holders.
otus = [] #list to fill with OTU id numbers
taxon_names = [] #list to fill with names of sequences that were matched to
otu_indices = [] #list to fill with tuples containing position to record
#each hit in the output table; row corresponds to otu and column corresponds to
#sample id
for i, line in enumerate(handle1) :
    #if read is hit or new seed
    if line[0] == "H" or line[0] =="S":
        if verbose :#print line number to stdout
            sys.stdout.write(str(i))
            sys.stdout.flush()
            sys.stdout.write("\b"*len(str(i)))
        separated = line.split("\t")
        if separated[3]>min_id:
            #get sample name
            id = separated[8].split(str(split_on))
            #build list of species' IDs
            if separated[1] not in otus :
                otus.append(separated[1])
            #build list of species' names
            #if hit, take last column as taxon name. split("/") used in case
            #matching record contains original sequence name, then the part
            #containing the read direction identifier and all after will be stripped   
            if line[0] == "H":
                matchcolumn=9
            #if new seed, take penultimate column as taxon name
            elif line[0] == "S":
                matchcolumn=8
            taxname=separated[matchcolumn].rstrip("\n")
            if taxname.split(str(split_on))[0] not in sample_ids:
                if db_name_part=="All":
                    pass
                else:
                    taxname=' '.join([taxname.split(str(db_name_split))[i] \
                    for i in db_name_part])
            if taxname not in taxon_names :
                    taxon_names.append(taxname)
            #create tuples of indices of OTU ID and corresponding sample name 
            otu_indices.append([taxon_names.index(taxname)+1, \
            sample_ids.index(id[0])])
        


if verbose:
    print "Generating OTU table from hit counts..."
#count all pairs of samples and OTUs
otu_table = zeros( (len(otus)+1, len(sample_ids)), dtype=object )
for index in otu_indices :
    otu_table[index[0], index[1]]+=1
#store sample names in table
for i, id_name in enumerate(sample_ids) :
    otu_table[0, i] = id_name
#store OTU ID in table
for i, otu in enumerate(otus) :
    otu_table[i+1, 0] = otu
#store OTU names
for i, tax_name in enumerate(taxon_names) :
    otu_table[i+1, -1] = tax_name
    

#Check whether a "Noise" row should be created
if noise:
    #store header
    temp_table=[otu_table[0]]
    noise_list = []
    for row in otu_table[1:] :
    #verify that sequence match was with database and not internal
        if max(row[1:-1]) > 0 and str(row[-1]).split(str(split_on))[0] not in sample_ids:
            temp_table.append(row)
    #otherwise, add to noise list
        if str(row[-1]).split(str(split_on))[0] in sample_ids:
             noise_list.append(row[1:-1])
    #clean up
    otu_table = None
    if len(noise_list) == 0:
        noise_table = zeros((1, len(sample_ids)-2), dtype=int)
    else:
        noise_table = array( noise_list , dtype=int )
    noise_list = None
    noise_id = array( ["0"], dtype=object)
    noise_name = array( ["Noise"], dtype=object)
    #number of noise reads in sample
    noise_counts = noise_table.sum(axis=0)
    noise_row = hstack( (noise_id, noise_counts, noise_name) )
    curated_otu_table=vstack( (temp_table, noise_row) )
    noise_table = None
    temp_table = None
else: #If no noise row required, use table currently stored in otu_table
    curated_otu_table = otu_table
if verbose:
    print "Writing data"
savetxt(handle2, curated_otu_table, fmt='%s', delimiter='\t')    
if verbose:
    print "Finished."
handle1.close()
handle2.close()
handle3.close()

