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
	#/usr/bin/time -o time.log -v python3 $testfile $lowdinbin > out.log  && echo -n $(cat out.log) ; grep "User time" time.log | gawk '{print "\t" $4 " Sec" }'
	#/usr/bin/time -o time.log -f "%e" python3 $testfile $lowdinbin > out.log  && echo -e -n $(cat out.log) ; echo -e $(cat time.log) " sec" | column -t -s$'\t'
	/usr/bin/time -o time.log -f "%e" python3 $testfile $lowdinbin > out.log  && printf "%-60s \t %s sec \n" "$(cat out.log)" $(cat time.log) 

    status=$((status + $?))
done
sh clean.sh
exit $status
