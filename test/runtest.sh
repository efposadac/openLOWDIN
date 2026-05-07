#!/bin/bash

if [ -z $1 ]; then
    EXENAME="openlowdin"
    if [ -e ../CONFIG ]; then
        EXENAME=`gawk '($1~/EXENAME/){print $3}' ../CONFIG`
    fi
else
    EXENAME=$1
fi

mkdir -p testResults_$EXENAME


date=$(date '+%Y-%m-%d_%H-%M-%S')
echo $date
echo "Testing with executable:" $EXENAME
echo "Saving outputs to " testResults_$EXENAME
RESULTS_LOG="testResults_$EXENAME/maketest_$date.log"

# copy fchk files. All tests/*fchk will be deleted with make clean
cp fchk/*fchk .

for testfile in `ls *.py`; do
    #Run test
    testName=`echo $testfile | gawk '{print substr($1,1,length($1)-3)}'`

    output=$( time -f "%e" -o time.log python3 "$testName.py" "$EXENAME" 2> error.log )
    error=$( cat error.log ) 
    duration=$( cat time.log ) 
 
    if [ -z "$error" ]; then
         printf "%-60s \t %s sec \n" "$output" "$duration" | tee -a "$RESULTS_LOG" 
    else 
         printf "%-60s \t %s sec \n %s \n" "$output" "$duration" "$error"  | tee -a "$RESULTS_LOG" 
    fi

    rm -f error.log 
    rm -f time.log 

    #Save results 
    find . -maxdepth 1 -name $testName.out -exec mv -t testResults_$EXENAME {} \;
    find . -maxdepth 1 -name $testName"*molden" -exec mv -t testResults_$EXENAME {} \;
    find . -maxdepth 1 -name $testName"*cub" -exec mv -t testResults_$EXENAME {} \;
    find . -maxdepth 1 -name $testName"*dens" -exec mv -t testResults_$EXENAME {} \;
    find . -maxdepth 1 -name $testName"*orb*" -exec mv -t testResults_$EXENAME {} \;
done

status=`grep -c "NOT OK" testResults_$EXENAME/maketest_$date.log`

if [ $status -gt 0 ]; then     
    echo $status "tests failed"
else
    echo "All tests completed successfully"
fi

exit $status
