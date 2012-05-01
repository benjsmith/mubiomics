#!/usr/bin/env python

#name2tax.py

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
from numpy import array, vstack, hstack, arange
from time import strftime



parser = argparse.ArgumentParser(description='''Map sample names in a csv file
to their taxonomic assignment at a chosen taxonomic rank. For use with
'taxtastic' reference packages.''')
parser.add_argument('-i','--infile', required=True, nargs=1,
                    type=str, help='''input filepath. Should be a tab-delimited
                    OTU table.''')
parser.add_argument('-o', '--outfile', required=True, nargs=1, type=str,
                    help='''output filepath. The resulting pooled OTU table
                    containing OTUs at desired taxonomic level is written here.
                    ''')
parser.add_argument('-t', '--taxfile', required=True, nargs=1, type=str,
                    help='''filepath to taxonomy table. Must be a csv format
                    table containing in the format created by Erick Matsen
                    et al's 'taxtastic'.
                    ''')
parser.add_argument('-m', '--mapfile', required=True, nargs=1, type=str,
                    help='''mapping file filepath. Must be a csv file
                    mapping from sequence name to NCBI taxon id, as created by
                    Erick Matsen et al's 'taxtastic'.
                    ''')
parser.add_argument('-r', '--rank', required=True, nargs=1, type=str,
                    help='''desired rank to map sequence names in OTU table to.
                    Must correspond to one of the rank names in the taxonomy
                    table, i.e., root, below_root, superkingdom,
                    below_superkingdom, below_below_superkingdom, superphylum,
                    phylum, subphylum, class, subclass, order, below_order,
                    below_below_order, suborder, family, below_family, genus,
                    species_group, species_subgroup, species, below_species.
                    Depending on available data in the taxonomy table, some
                    ranks may not have a mapping.''')

args = parser.parse_args()
infile = args.infile[0]
outfile = args.outfile[0]
taxfile = args.taxfile[0]
mapfile = args.mapfile[0]
rank = args.rank[0] 



print "\nRun started " + strftime("%Y-%m-%d %H:%M:%S") + "."

inhandle = csv.reader(open(infile, 'rU'), delimiter='\t')
outhandle = csv.writer(open(outfile, 'wb'), delimiter='\t')
#Write header lines of output file.
outhandle.writerow([''.join(inhandle.next()) + ', Rank: ' + rank])
outhandle.writerow(inhandle.next())

#collect sample names, using first line of file
inhandle = csv.reader(open(infile, 'rU'), delimiter='\t')
inhandle.next()
sample_ids = [column for column in \
inhandle.next() if \
re.search(column, "#OTU ID"'|'"Consensus Lineage")==None]

#make matrix from input file
inhandle = csv.reader(open(infile, 'rU'), delimiter='\t')
in_mat = [line for i, line in enumerate(inhandle) if i>1]


#make mapping dictionary from sequence name to taxonomic id
maphandle = csv.reader(open(mapfile, 'rU'))
map_dict = dict([(line[0], line[2]) \
    for i, line in enumerate(maphandle) if i>0])
#make mapping dictionary from tax_id to tax_name
taxhandle = csv.reader(open(taxfile, 'rU'))
tax_dict = dict([(line[0], line[3]) \
    for i, line in enumerate(taxhandle) if i>0])
#make matrix from taxonomy file for lookup
taxhandle = csv.reader(open(taxfile, 'rU'))
tax_mat = [line for line in taxhandle]

for pos in arange(len(in_mat)):
    if in_mat[pos][-1] != "Noise":
        row = [part[0] for part in tax_mat].index(map_dict[in_mat[pos][-1]])
        col = tax_mat[0].index(rank)
        if tax_mat[row][col] == '':
            in_mat[pos][-1] = \
            str(tax_dict[map_dict[in_mat[pos][-1]]]) +": no "+rank+" data"
        else:
            in_mat[pos][-1] = tax_dict[tax_mat[row][col]]
    else:
        pass
for row in in_mat:
    outhandle.writerow(row)

print "Run finished " + strftime("%Y-%m-%d %H:%M:%S") + "."
