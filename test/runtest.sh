#!/bin/bash/

for testfile in `ls *.py`; do
	python $testfile
    status=$((status + $?))
done

exit $status