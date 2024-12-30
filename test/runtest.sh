#!/bin/bash

if [ -z $1 ]
then
    lowdinbin="lowdin2"
else
    lowdinbin=$1
fi
echo $lowdinbin

for testfile in `ls *.py`; do
    #echo $testfile
	#python3 $testfile $lowdinbin
	#/usr/bin/time -v python3 $testfile $lowdinbin 2>&1 | (echo -n $(head -1) ; grep "User time" | gawk '{print "\t" $4 " Sec" }' )
	/usr/bin/time -o time.log -v python3 $testfile $lowdinbin > out.log  && echo -n $(cat out.log) ; grep "User time" time.log | gawk '{print "\t" $4 " Sec" }'

    status=$((status + $?))
done
sh clean.sh
exit $status
