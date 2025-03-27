#!/bin/bash

if [ -z $1 ]; then
    lowdinbin="openlowdin"
    if [ -e ../CONFIG ]; then
        lowdinbin=`gawk '($1~/EXENAME/){print $3}' ../CONFIG`
    fi
else
    lowdinbin=$1
fi
echo "Testing with lowdinbin:" $lowdinbin

for testfile in `ls *.py`; do
    python3 $testfile $lowdinbin
    status=$((status + $?))
done

if [ $status -gt 0 ]; then     
    echo $status "tests failed"
else
    echo "All tests completed successfully"
fi

exit $status
