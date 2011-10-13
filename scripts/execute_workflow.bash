#!/bin/bash
#$ -S /bin/bash -cwd -j y -l mf=2G -N demult_execute_worfklow -V
SPF=$(( ${SGE_TASK_ID}-1 ))
STN=$(( ${SPF}*10000000 ))
if [ "$SGE_TASK_ID" -le "10" ]; then
	SPF="0$SPF"
fi
python ${SCRIPTLOC}/demultiplexer.py <(echo -e "-i${1}.split${SPF}\n-o${1}.out${SPF}\n-n$STN" | cat - $2)
${SCRIPTLOC}/usearch --sort ${1}.split$SPF --output ${1}.sorted$SPF
${SCRIPTLOC}/usearch --cluster ${1}.sorted$SPF --db $3 --uc ${1}.uc$SPF --id 0.89 --rev --iddef 2 --local
python ${SCRIPTLOC}/ucstripper.py ${1}.uc$SPF ${1}.table$SPF $4
