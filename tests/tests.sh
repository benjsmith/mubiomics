#!/bin/bash
#HPV genotyping workflow as test shell script. 

#Data was generated from an Illumina HiSeq2000 sequencing run of multiple
#clinical samples containing HPV infections of different types.

#Run qc on fastq file using the command line arguments stored in
#qc_args.txt.

qc.py @qc_args.txt

#Run demultiplexing on QCed fastq file using the command line arguments stored
#in dm_args.txt.

demultiplexer.py @dm_args.txt

#Run usearch on assigned reads against HPV database. usearch program must be
#on shell $PATH for this to work.
#To download usearch visit http://www.drive5.com/usearch/

#First sort sequences. 
usearch --sort testout.fasta --output testout.sorted.fasta
#Then run search against a database on sorted sequences.
usearch --cluster testout.sorted.fasta --db PV.MYn.db.fasta --uc \
testout.uc --id 0.89 --rev --iddef 2 --global

#Run uclust output file parser to create table of results.
ucstripper.py -i testout.uc -o table.txt -m testmap.txt

#Optionally, pool OTUs with identical names and remove OTUs with less than a 
#minimum proprotion of counts in any sample.
pool_otus.py table.txt pooled_table.txt -m 0.01 -r

