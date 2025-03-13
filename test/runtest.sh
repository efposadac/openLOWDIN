#!/bin/bash

if [ -z $1 ]
then
    lowdinbin="openlowdin"
else
    lowdinbin=$1
fi
echo $lowdinbin

for testfile in `ls *.py`; do
    #echo $testfile
    python3 $testfile $lowdinbin
    status=$((status + $?))
done

exit $status
