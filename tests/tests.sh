#!/bin/bash
echo ""
echo "#######################################"
echo "Sample workflow as test shell script."
echo "#######################################"
echo ""
echo "These operations should create a new directory called qc and"
echo "run to completion without errors, if not, please email me."
echo "benjamincharlessmith@gmail.com"
echo ""
echo "150,000 reads are provided in the test data and will take about 25 mins"
echo "to complete; the perfect opportunity for a break!"
echo "If you don't want to wait that long, you could save a subset of the reads"
echo "and replace testdata.fastq with the subset."
echo ""
echo "Run qc on test fastq file using the command line arguments stored in"
echo "qc_args.txt."
echo ""
echo ">qc.py @qc_args.txt"
qc.py @qc_args.txt
echo""
echo "Run demultiplexing on QCed fastq file using the command line arguments"
echo "stored in dm_args.txt."
echo ""
echo ">demultiplexer.py @dm_args.txt"
demultiplexer.py @dm_args.txt
echo ""
echo ""
echo "Run usearch on assigned reads against HPV database. usearch program must"
echo "be on shell \$PATH for this to work."
echo "To download usearch visit http://www.drive5.com/usearch/"
echo ""
echo "First, sort sequences."
echo ""
echo ">usearch --sort qc/testout.fasta --output qc/testout.sorted.fasta"
usearch --sort qc/testout.dmx.fasta --output qc/testout.sorted.fasta
echo "Then, run usearch against a database on the sorted sequences."
echo ""
echo ">usearch --cluster qc/testout.sorted.fasta --db testdb.fasta --uc \
qc/testout.uc --id 0.89 --rev --iddef 2 --global"
usearch --cluster qc/testout.sorted.fasta --db testdb.fasta --uc \
qc/testout.uc --id 0.89 --rev --iddef 2 --global
wait
echo ""
echo "Run uclust output file parser to create table of results."
echo ""
echo ">ucstripper.py -i qc/testout.uc -o qc/table.txt -m testmap.txt -v"
ucstripper.py -i qc/testout.uc -o qc/table.txt -m testmap.txt -v
wait
echo ""
echo "Optionally, pool OTUs with identical names and remove OTUs with less"
echo "than a minimum proprotion of counts in any sample."
echo ""
echo ">pool_otus.py qc/table.txt qc/pooled_table.txt -m 0.01 -r"
pool_otus.py -i qc/table.txt -o qc/pooled_table.txt -m 0.01 -r
echo ""
echo "Finished."

