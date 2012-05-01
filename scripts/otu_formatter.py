#!/usr/bin/env python
#otu_formatter.py

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
from numpy import zeros, arange, hstack, vstack, array
from itertools import product
from time import strftime
import cPickle as pickle


parser = argparse.ArgumentParser(description='''From multiple input OTU tables,
                    produce the same number of output OTU tables where each
                    new table contains the same number of rows and columns,
                    the same row and column headings and the same ordering along
                    rows and columns.''')
parser.add_argument('-i','--infiles', required=True, nargs='+',
                    type=str, help='''input file list. Should be one or more
                    tab-delimited OTU tables.''')

#Setting required variables from command line arguments.
args = parser.parse_args()
infiles = args.infiles

print "\nRun started " + strftime("%Y-%m-%d %H:%M:%S") + "."

all_ids = []#list to contain all observed sample ids
all_otus = []#list to contain all observed taxon names
all_dicts = []#list to contain a dictionary for each file, mapping taxon names
# and samples to counts
outhandles=[]#initialise list that will contain csv.writer objects
for i, file in enumerate(infiles):
    file_dict = {}#dictionary to map taxon and sample names to a count
    inhandle = csv.reader(open(file, "rU"), delimiter='\t')
    outfile = str(os.path.splitext(file)[0]) + "_final" + \
    str(os.path.splitext(file)[1])#get outpath from inpath
    outhandles.append(csv.writer(open(outfile, "wb"), delimiter='\t'))# output file handle
    #write first header line of output file
    outhandles[i].writerow([''.join(inhandle.next()) + ', formatted'])
    #get sample ids in current file from second line in file
    idline = inhandle.next()
    sample_ids = [column for column in idline if \
re.search(column, "#OTU ID"'|'"Consensus Lineage")==None]
    #scan through current file ids and add to all_ids any that are not already contained.
    for id in sample_ids:
        if id not in all_ids:
            all_ids.append(id)
    #scan through each row in the remaining input file, which is the matrix of counts
    for row in inhandle:
        #if the taxon name is not already present, add to all_otus
        if row[-1] not in all_otus:
            all_otus.append(row[-1])
        #scan through each cell containing counts and add data to file_dict,
        #where the key is a pickled list, [taxon name, sample id]
        # corresponding to that cell and the value is the count in that cell
        for i, cell in enumerate(row[1:-1]):
            file_dict[pickle.dumps([row[-1], sample_ids[i]])] = cell 
    #append the completed file_dict to the all_dicts
    all_dicts.append(file_dict)

#sort taxon name and sample id lists
all_ids.sort()
all_otus.sort()
#create template matrix of zeros with dimensions taken from all_ids and all_otus.
temp_mat = zeros((len(all_otus), len(all_ids)), dtype=int)

for i, file in enumerate(infiles):
    outhandle = outhandles[i]
    file_mat = temp_mat.copy()
    for key, value in all_dicts[i].items():
        row = all_otus.index(pickle.loads(key)[0])
        col = all_ids.index(pickle.loads(key)[1])
        file_mat[row, col] = value
    idline = ["#OTU ID"] + all_ids + ["Consensus Lineage"]
    out_mat = hstack((array([[i]for i in arange(len(all_otus))]), \
                file_mat, array([[name] for name in all_otus])))
    final_mat = vstack((idline, out_mat))
    for row in final_mat:
        outhandle.writerow(row)

print "Run finished " + strftime("%Y-%m-%d %H:%M:%S") + "."
