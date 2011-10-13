#!/bin/bash
set -o errexit -o nounset -o pipefail -o allexport
PATH=/share/apps/python/current/bin:$PATH
if [ "$#" -ne "4" ]; then
	echo "usage: demultiplexer.bash <input.fastq> <args.txt> <db.txt> <map.txt>" >&2
	exit 1
fi
SCRIPTLOC=$(cd "`dirname \"$0\"`"; echo "$PWD")
qsub -sync y split.bash $1
qsub -sync y -t 1-`ls -1 ${1}.split* | wc -l` execute_workflow.bash $@
qsub -sync y merge.bash $1
