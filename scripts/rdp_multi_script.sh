#!/usr/bin/env bash
#Script to run pplacer on multiple sample files.
#Usage: rdp_multi_script.sh

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

SAVEIFS=$IFS

IFS=$(echo -en "\n\b")

echo "###########################"
echo "##  rdp_multi_script.sh  ##"
echo "###########################"
echo ""
echo "This script will run rdp_multiclassifier.jar for a fasta file containing \
multiple sequences from multiple samples, against a custom trained library."
echo "It will prompt you for the following inputs:"
echo "-the location of the rRNAClassifier.properties file for you library (the enclosing \
folder should contain all additional library files produced when training your library)."
echo "-the path to the file you wish to save assignment results to (should end in .csv)."
echo "-the path to the fasta file containing your sequences."
echo "-the minimum bootstrap confidence at which you wish to accept classifications (for reads below \
250 bp, it is suggested that you use a value of 0.5)"
echo ""
echo "##########################"
echo ""


read -p "Enter path to rRNAClassifier.properties file and press \
return: " PROPFILE
read -p "Enter path to file you wish to save results in and press \
return: " OUTFILE
read -p "Enter path to fasta file containing sequences for classification and press \
return: " SEQFILE
read -p "Enter minimum bootstrap confidence: " CONF
echo "##########################"
echo ""

PROPFILE=$(echo "${PROPFILE}" | sed -e 's/^ *//g;s/ *$//g')
OUTFILE=$(echo "${OUTFILE}" | sed -e 's/^ *//g;s/ *$//g')
SEQFILE=$(echo "${SEQFILE}" | sed -e 's/^ *//g;s/ *$//g')


echo ""
echo "Running rdp_multiclassifier on input file."
echo "Press ^C (Control+C) to prematurely terminate the run at any time."
echo ""
java -Xmx1g -jar /Applications/rdp_multiclassifier/dist/multiclassifier.jar \
--train_propfile="${PROPFILE}" --conf="${CONF}" --format=fixRank \
--assign_outfile="${OUTFILE}" "${SEQFILE}"
echo ""
echo ""
echo "Run finished."
echo "Results can be found in "${OUTFILE}""

IFS=$SAVEIFS