#!/usr/bin/env bash
#Shell script for running sqlite3 commands.
#Usage: sqlite_script.sh

read -p "Enter path to database file: " DB
read -p "Enter minimum acceptable likelihood (0-1): " THRESHOLD
read -p "Enter taxonomic rank at which to extract records: " TAX

sqlite3 -header -csv "${DB}" "
SELECT pc.placement_id,
       placement_names.origin,
       taxa.tax_name,
       pc.desired_rank,
       pc.likelihood
FROM placement_classifications AS pc
INNER JOIN taxa
ON pc.tax_id=taxa.tax_id
INNER JOIN placement_names
ON pc.placement_id=placement_names.placement_id
WHERE pc.likelihood>"${THRESHOLD}"
AND pc.desired_rank='"${TAX}"'
AND pc.rank='"${TAX}"'
ORDER BY pc.placement_id    
" >"${TAX}"_at_"${THRESHOLD}".csv