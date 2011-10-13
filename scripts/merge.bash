#!/bin/bash
#$ -S /bin/bash -cwd -j y -l mf=2G -N demult_merge -V
set -o errexit -o nounset
python ${SCRIPTLOC}/pool_otus.py <(cat ${1}.table* | grep -v '^#') ${1}.pooled_table.txt -m 5 -r
