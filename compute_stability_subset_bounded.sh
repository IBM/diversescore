#!/bin/bash

# $1 domain
# $2 problem
# $3 plans folder
# $4 number of plans in the folder (N)
# $5 number of plans to choose (M<=N)
# $6 bound for the metric
# $7 (optional) plan file name - default is sas_plan
## Assuming plan files are named <plan file name>.1 ... <plan file name>.N 

if [ "$#" -lt 6 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

if [ "$#" -gt 7 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

if [ "$5" -gt "$4" ]; then
    echo "Illegal parameter - subset size cannot be greater than the set size"
    exit 1
fi

FNAME=""
if [ "$#" -eq 7 ]; then
    FNAME="--internal-plan-file "$7
fi

SCORE="subset_bounded(compute_stability_metric=true,aggregator_metric=min,plans_as_multisets=true,plans_subset_size=$5,metric_bound=$6,metric_computation_method=mip,dump_plans=true)"

SOURCE="$( dirname "${BASH_SOURCE[0]}" )"
$SOURCE/fast-downward.py $1 $2 --diversity-score $SCORE --internal-plan-files-path $3 --internal-num-plans-to-read $4 $FNAME


