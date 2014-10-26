#!/bin/bash/

for testfile in `ls *.py`
do
	python $testfile
done
