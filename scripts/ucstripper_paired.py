#!/usr/bin/env python
#ucstripper_paired.py
#convert results from a usearch database query to an OTU table for paired
#end data.

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
    description='''Strip data from a paired-end usearch .uc output file and
    produce an OTU table. The usearch input file must have been
    demultiplexed.''',
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
                    headings in the output. Should contain both forward and
                    reverse names. If sample names in the mapping file
                    end with '.F' or '.R', respectively, set the -s or
                    --split_on parameter to a period.''')
parser.add_argument('-s', '--split_on', nargs=1, default=['.'], type=str,
                    help='''set the chracter to split the sample name on. E.g.,
                    if the read name is sample1_00001 and the split is set to
                    "_", the sample name will be read as "sample1".''')
parser.add_argument('-n', '--noise', action='store_true',
                    help='''condense counts for all clusters that do not match
                    any database sequences into a group called "Noise".''')
parser.add_argument('-p', '--paired_only', action='store_true',
                    help='''only count taxa where reads from both directions
                    agree. If left unset, reads without partner hits will also
                    be counted.''')
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
verbose = args.verbose
noise = args.noise
paired_only = args.paired_only

#open input/output handles
handle1 = open(infile, "rb") # uclust output file handle
handle2 = open(outfile, "wb") # output text file handle
handle3 = open(mapfile, 'rb') #mapping file handle

print "Running..."

handle2.write("#ucstripper OTU Table\n") # add header line to outout file

#parse mapping file to build list of sample names, these will be the column
#headers for the count data.
sample_ids = ['#OTU ID']
for line in handle3 :
    if line[0] != "#" :
        sample_ids.append(line.split()[0].split(split_on)[0])
sample_ids.append("Consensus Lineage")

if verbose:
    print "Processing reads for " + str(len(sample_ids)-2) + " samples."
    for id in sample_ids:
        print id

#parse the input, 'hit' by 'hit', and add info from columns to
#data holders.
otus = []
taxon_names = []
otu_indices = []
reads = {}
read2sample = {}
for i, line in enumerate(handle1) :
    #if hit
    if line[0] == "H" or line[0] =="S" :
        if verbose :#print line number to stdout
            sys.stdout.write(str(i))
            sys.stdout.flush()
            sys.stdout.write("\b"*len(str(i)))
        separated = line.split("\t")
        #compare original read name to processed_names to determine wheher to
        #process this line
        try:
            orig = str(separated[8].split()[1].split('/')[0]) #look for
            #original name at second position in column 8 and strip e.g,
            #'/2' from end.
        except IndexError:
            print "Error produced by read names. Demultiplexed read names \
should contain the original read name after the new sample name, \
separated by a space."
            sys.exit(1)
        #get sample name
        id = separated[8].split(str(split_on))
        #build list of species' IDs
        if separated[1] not in otus :
            otus.append(separated[1])
        #build dict mapping read name to sample_ids index. Used for flushing
        #reads dictionary
        read2sample[orig] = sample_ids.index(id[0])
        #build list of species' names
        #if hit, take last column as taxon name. split("/") used in case
        #matching record contains original sequence name, then the part
        #containing the read direction identifier and all after will be stripped
        if line[0] == "H":
            tax_name = separated[9].rstrip("\n").split("/")[0]
        #if new seed, take penultimate column as taxon name
        elif line[0] == "S":
            tax_name = separated[8].split("/")[0]
        if tax_name not in taxon_names :
            taxon_names.append(tax_name)
        try:
            #build dictionary containing key = original read name, value =
            #taxon it is associated with.
            if reads[orig] == tax_name:
            # if existing taxon name is the same as the current one, create
            #tuple containing indices which will be used to fill otu_table
            #with count data by position.
                otu_indices.append([otus.index(separated[1])+1, \
                sample_ids.index(id[0])])
                #also delete this dictionary entry, as the pair is found and
                #agrees in the taxon assignment
                del reads[orig]
            else:#if there is no agreement,
                #just delete the dictionary entry, as the reads should not
                #be counted
                del reads[orig]
        except KeyError:
            reads[orig] = tax_name
            

if verbose:
    print "All matching pairs have been counted."
    print "Flushing list of unpaired reads..."

#parse the remaining reads, these will be singletons, and add their otu
#indices to the list, only if -p not specified.
if not paired_only:
    for k, v in reads.iteritems():
        otu_indices.append([taxon_names.index(v)+1, \
        read2sample[k]])#+1 because of header line in table
    reads={}
if verbose:
    print "Generating OTU table from hit counts..."

#count all pairs of samples and OTUs based on positional tuples in otu_indices
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
#store header
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
savetxt(handle2, curated_otu_table, fmt='%s', delimiter='\t')    
print "Finished."
handle1.close()
handle2.close()
handle3.close()

