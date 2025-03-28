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

for testfile in `ls *.py`; do
    #Run test
    testName=`echo $testfile | gawk '{print substr($1,1,length($1)-3)}'`
    python3 $testName.py $EXENAME | tee -a testResults_$EXENAME/maketest_$date.log
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
