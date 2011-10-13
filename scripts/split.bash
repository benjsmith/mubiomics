#!/bin/bash
#$ -S /bin/bash -cwd -j y -N demult_split -l mf=2G -V
set -o errexit -o nounset
split -dl10000000 $1 ${1}.split
