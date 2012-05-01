#!/usr/bin/env bash
#Script to run pplacer on multiple sample files.
#Usage: pplacer_batch.sh

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
echo "##   pplacer_batch.sh    ##"
echo "###########################"
echo ""
echo "This script will run a batch phylogenetic pipeline on a folder of fasta \
alignment files. The pipeline consists of the following steps:"
echo "-Deduplicate input sequences (i.e., remove any duplicate sequences)."
echo "-Run pplacer to determine phylogenetic placement of reads."
echo "-Reduplicate placement files (i.e., repopulate with duplicate information)."
echo "-Produce phylogentic trees for visualisation of sample composistion \
using guppy fat."
echo ""
echo "##########################"
echo ""


read -p "Enter directory containing input files and press \
return (no trailing slash): " INDIR
read -p "Enter directory containing reference package and press \
return (no trailing slash): " REFPKG

echo "##########################"
echo ""


INDIR=$(echo "${INDIR}" | sed -e 's/^ *//g;s/ *$//g')
REFPKG=$(echo "${REFPKG}" | sed -e 's/^ *//g;s/ *$//g')

DEDUPDIR=${INDIR%/*}/dedup_files

if [ ! -d "${DEDUPDIR}" ];
    then
        mkdir "${DEDUPDIR}"
        echo "dedup output directory created: $DEDUPDIR"
        echo ""
fi

OUTDIR=${INDIR%/*}/pplacer_out

if [ ! -d "${OUTDIR}" ];
    then
        mkdir "${OUTDIR}"
        echo "pplacer output directory created: $OUTDIR"
        echo ""
fi

GPYOUTDIR=${INDIR%/*}/guppy_fat_out

if [ ! -d "${GPYOUTDIR}" ];
    then
        mkdir "${GPYOUTDIR}"
        echo "guppy output directory created: $GPYOUTDIR"
        echo ""
fi


FILES=($(ls "${INDIR}"/*.fasta))

echo ""
echo "Running pplacer on input files."
echo "Press ^C (Control+C) to prematurely terminate the run at any time."
echo ""
for INFILE in ${FILES[@]}; do
    FILEID=${INFILE##*/}
    if [ ! -f "${OUTDIR}"/${FILEID%.*}.json ]; then
        echo "Processing ${FILEID}"
        seqmagick convert --deduplicate-sequences \
        --deduplicated-sequences-file="${DEDUPDIR}"/${FILEID%.*}.dedup \
        "${INFILE}" "${DEDUPDIR}"/${FILEID%.*}.deduped.fasta
        wait
        pplacer -c "${REFPKG}" --out-dir "${DEDUPDIR}"  "${DEDUPDIR}"/${FILEID%.*}.deduped.fasta
        wait
        guppy redup -d "${DEDUPDIR}"/${FILEID%.*}.dedup "${DEDUPDIR}"/${FILEID%.*}.deduped.jplace \
        -o "${OUTDIR}"/${FILEID%.*}.json
        wait
        rm "${DEDUPDIR}"/${FILEID%.*}.deduped.jplace
        wait
        guppy fat -c "${REFPKG}" "${OUTDIR}"/${FILEID%.*}.json -o "${GPYOUTDIR}"/${FILEID%.*}.xml
        wait
        echo ""
    fi
    done

echo "Batch run finished."
echo "Results can be found in:"
echo "${DEDUPDIR}"
echo "${OUTDIR}"
echo "${GPYOUTDIR}"

IFS=$SAVEIFS