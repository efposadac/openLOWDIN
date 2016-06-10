#!/bin/bash

for testfile in `ls *.py`; do
    echo $testfile
	python $testfile
    status=$((status + $?))
done
sh clean.sh
exit $status