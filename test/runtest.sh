#!/bin/bash/

for testfile in `ls *.py`; do
	python $testfile
    status=$((status + $?))
done
sh clean.sh
exit $status