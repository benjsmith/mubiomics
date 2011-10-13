#!/usr/bin/env python

#Author: Benjamin Smith, 2011

from Bio import SeqIO
from Bio.Seq import Seq
import sys, os, re, argparse, csv
from numpy import array, vstack, hstack, arange
from time import strftime



parser = argparse.ArgumentParser(description='''Map sample names in the header
lines of a fasta file to their taxonomic assignment at a chosen taxonomic rank.
For use with 'taxtastic' reference packages.''')
parser.add_argument('-i','--infile', required=True, nargs=1,
                    type=str, help='''input filepath. Should be a FASTA format
                    file.''')
parser.add_argument('-o', '--outfile', required=True, nargs=1, type=str,
                    help='''output filepath. The resulting FASTA file with
                    taxonomic names in the headers is written here.
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

inhandle = open(infile, 'rU')
record_iter = SeqIO.parse(inhandle, "fasta")
try :
    outhandle = open(outfile, "wb")
except IOError:
    outdir = str(os.path.dirname(os.path.abspath(outfile)))
    os.mkdir(outdir)
    outhandle = open(outfile, "wb") # create output fasta/fastq file if it
    #doesn't exist


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

for seq in record_iter:
    new_seq = seq
    new_seq.name = ""
    new_seq.description = ""
    row = [part[0] for part in tax_mat].index(map_dict[seq.id])
    col = tax_mat[0].index(rank)
    if tax_mat[row][col] == '':
        new_seq.id = seq.id + "|" + \
        str(tax_dict[map_dict[seq.id]]) +": no "+rank+" data"
    else:
        new_seq.id = seq.id + "|" + tax_dict[tax_mat[row][col]]
    #print new_seq
    SeqIO.write(new_seq, outhandle, "fasta")

#for pos in arange(len(in_mat)):
#    if in_mat[pos][-1] != "Noise":
#        row = [part[0] for part in tax_mat].index(map_dict[in_mat[pos][-1]])
#        col = tax_mat[0].index(rank)
#        if tax_mat[row][col] == '':
#            in_mat[pos][-1] = \
#            str(tax_dict[map_dict[in_mat[pos][-1]]]) +": no "+rank+" data"
#        else:
#            in_mat[pos][-1] = tax_dict[tax_mat[row][col]]
#    else:
#        pass
#for row in in_mat:
#    outhandle.writerow(row)

print "Run finished " + strftime("%Y-%m-%d %H:%M:%S") + "."
