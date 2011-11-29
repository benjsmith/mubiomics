#!/usr/bin/env bash
#Shell script for running sqlite3 commands.
#Usage: sqlite_script.sh

read -p "Enter path to database file: " DB

OUTPATH=$(dirname "${DB}")

sqlite3 -header -csv "${DB}" "
SELECT placement_names.name,
       taxa.tax_name,
       pc.rank,
       pc.likelihood
FROM placement_classifications AS pc
INNER JOIN taxa
ON pc.tax_id=taxa.tax_id
INNER JOIN placement_names
ON pc.placement_id=placement_names.placement_id
WHERE pc.desired_rank=pc.rank
ORDER BY placement_names.name    
" >"${OUTPATH}"/raw_results_table.csv