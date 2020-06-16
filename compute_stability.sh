#!/bin/bash

# $1 domain
# $2 problem
# $3 plans folder
# $4 number of plans in the folder (N)
# $5 (optional) plan file name - default is sas_plan
## Assuming plan files are named <plan file name>.1 ... <plan file name>.N 

if [ "$#" -lt 4 ]; then
    echo "Illegal number of parameters"
fi

if [ "$#" -gt 5 ]; then
    echo "Illegal number of parameters"
fi

FNAME=""
if [ "$#" -eq 5 ]; then
    FNAME="--internal-plan-file "$5
fi

SOURCE="$( dirname "${BASH_SOURCE[0]}" )"
$SOURCE/fast-downward.py $1 $2 --diversity-score "score(compute_stability_metric=true,aggregator_metric=avg,plans_as_multisets=false)" --internal-plan-files-path $3 --internal-num-plans-to-read $4 $FNAME
