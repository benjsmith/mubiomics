#!/usr/bin/env python
#taxtastic2rdp.py
#converts a taxtastic database to one suitable for use with RDP Classifier (Wang et al., 2010)

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
from Bio import SeqIO



parser = argparse.ArgumentParser(description='''Produce an RDP
(multi-)classifier reference package from a taxtastic reference package.''')
parser.add_argument('-f','--fasta_in', required=True, nargs=1,
                    type=str, help='''input filepath. Should be a fasta file
                    containing the sequences to convert to an RDP classifier
                    reference package.''')
parser.add_argument('-t', '--tax_in', required=True, nargs=1, type=str,
                    help='''filepath to taxonomy table. Must be a csv format
                    table containing in the format created by Erick Matsen
                    et al's 'taxtastic'.
                    ''')
parser.add_argument('-m', '--map_in', required=True, nargs=1, type=str,
                    help='''mapping file filepath. Must be a csv file
                    mapping from sequence name to NCBI taxon id, as created by
                    Erick Matsen et al's 'taxtastic'.
                    ''')
parser.add_argument('-F', '--fasta_out', required=True, nargs=1, type=str,
                    help='''fasta library output filepath. The resulting RDP
                    format fasta reference file is written here.
                    ''')
parser.add_argument('-T', '--tax_out', required=True, nargs=1, type=str,
                    help='''taxonomy output filepath. The resulting RDP format
                    taxonomy file is written here.
                    ''')
parser.add_argument('-d', '--depth', nargs=1, default=['genus'], type=str,
                    help='''maximum desired depth.
                    ''')


args = parser.parse_args()
fasta_in = args.fasta_in[0]
tax_in = args.tax_in[0]
map_in = args.map_in[0]
fasta_out = args.fasta_out[0]
tax_out = args.tax_out[0]
max_depth = args.depth[0]



print "\nRun started " + strftime("%Y-%m-%d %H:%M:%S") + "."

####################################
#Create RDP taxonomy file
###################################

#make dictionary mapping rank to phylogentic depth as integer e.g., {root: 0}
tax_inhandle = csv.reader(open(tax_in, 'rU'))
taxdepth = {}
ranklist = []
for i, rank in enumerate(tax_inhandle.next()[4:]):
    taxdepth[rank] = i
    ranklist.append(rank)
depths = sorted(taxdepth.items(), key=lambda (k,v): v)
#make list of all tax_ids so that can re-number any that are not integers, 
#without conlicting with existing tax_ids
depth_num = taxdepth[max_depth]
all_tax_ids = [int(line[0].split("_")[0]) for line in tax_inhandle]
all_tax_ids.sort()
#create dictionary containing mappings betwen old and new tax_ids
tax_id_map = {}
#make rdp taxonomy list from taxtastic taxonomy file
tax_inhandle = csv.reader(open(tax_in, 'rU'))#reopen tax_inhandle as it was
#exhausted by the all_tax_ids list comprehension.
tax_inhandle.next() #skip header line
rdp_taxlist = []
for i, line in enumerate(tax_inhandle):
    if taxdepth[line[2]] <= depth_num:#select only lines at "genus" level or higher
        rel_depth = len([part for part in line[5:] if part != ''])
        if len(line[0].split("_")) > 1:#select ids containing '_' for tax_id reassignment
            taxid = all_tax_ids[-1] + 1 #new tax_id is 1 + the largest observed one
            tax_id_map[line[0]] = taxid#add new one to mapping dictionary
            all_tax_ids.append(taxid)
            rdp_taxlist.append([int(taxid), line[3], tax_id_map[line[1]], \
                                int(rel_depth), line[2]])
        else:
            tax_id_map[line[0]] = int(line[0])
            rdp_taxlist.append([int(line[0]), line[3], tax_id_map[line[1]], \
                                int(rel_depth), line[2]])

#write RDP taxonomy file
tax_outhandle = csv.writer(open(tax_out, 'w'), delimiter='*')
for row in rdp_taxlist:
    tax_outhandle.writerow(row)

##################################
#create fasta file of reference sequences with additional taxonomic
#information in the headers
#################################
fasta_inhandle = open(fasta_in, 'rU')
record_iter = SeqIO.parse(fasta_inhandle, 'fasta')

#open output file for writing
try :
    fasta_outhandle = open(fasta_out, "wb")
except IOError:
    outdir = str(os.path.dirname(os.path.abspath(fasta_out)))
    os.mkdir(outdir)
    fasta_outhandle = open(fasta_out, "wb") # output fasta file

#make mapping dictionary from sequence name to taxonomic id
maphandle = csv.reader(open(map_in, 'rU'))
map_dict = dict([(line[0], line[2]) \
    for i, line in enumerate(maphandle) if i>0])
#make mapping dictionary from tax_id to tax_name
taxhandle = csv.reader(open(tax_in, 'rU'))
tax_dict = dict([(line[0], line[3]) \
    for i, line in enumerate(taxhandle) if i>0])
#make mapping from tax_name to tax_id
taxhandle = csv.reader(open(tax_in, 'rU'))
name_dict = dict([(line[3], line[0]) \
    for i, line in enumerate(taxhandle) if i>0])
#make matrix from taxonomy file for lookup
taxhandle = csv.reader(open(tax_in, 'rU'))
tax_mat = [line for line in taxhandle]

for rec in record_iter:
    newrec = rec
    newdesc = "root"
    row = [part[0] for part in tax_mat].index(map_dict[rec.id])
    depth_ctr = 0
    rel_depth = 0
    for rank in ranklist[1:]:#iterate through ranks after 'root'
        depth_ctr += 1
        col = tax_mat[0].index(rank)
        if tax_mat[row][col] == '':
            if depth_ctr == depth_num:#If info at specified depth does not exist, create an artificial taxon.
                parent_name = newdesc.split(";")[-1]
                artificial_genus = "genus below " + parent_name
                newdesc = newdesc + ";" + artificial_genus
                candidate_tax = [artificial_genus, \
                               name_dict[parent_name], \
                                rel_depth+1, max_depth] 
                for tax_row in rdp_taxlist:#test that new artificail genus is not in taxa already written
                    if candidate_tax == tax_row[1:]:
                        tax_flag = 1
                        break
                    else:
                        tax_flag = 0
                if tax_flag == 0:
                    taxid = all_tax_ids[-1] + 1 #new tax_id is 1 + the largest observed one
                    tax_id_map[taxid] = taxid#add new one to mapping dictionary
                    all_tax_ids.append(taxid)
                    name_dict[artificial_genus] = taxid
                    tax_dict[taxid] = artificial_genus
                    new_tax_row = [taxid] + candidate_tax
                    rdp_taxlist.append(new_tax_row)
                    tax_outhandle.writerow(new_tax_row)
                break
            else:
                pass
            #newdesc = newdesc + ";" + newdesc.split(";")[-1]
        else:
            rel_depth += 1
            newdesc = newdesc + ";" + tax_dict[tax_mat[row][col]]
            if depth_ctr == depth_num:
                break
            else:
                pass
    newrec.description = newdesc
    SeqIO.write(newrec, fasta_outhandle, 'fasta')
    

print "\nRun finished " + strftime("%Y-%m-%d %H:%M:%S") + ".\n"
print "Run the following command to train RDP classifier for your database:"
print "java -Xmx1g -cp rdp_classifier-2.3.jar \
edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker" \
+ tax_out + fasta_out + " 1 v1 test" + \
str(os.path.dirname(os.path.abspath(fasta_out)))
