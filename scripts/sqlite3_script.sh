#!/usr/bin/env bash
#Shell script for running sqlite3 commands.
#Usage: sqlite3_script.sh

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