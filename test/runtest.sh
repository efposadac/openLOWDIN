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
echo "| ----------------------------------- | ------- | --------------------------------------- |  "
echo "| Testname                            | Time(s) | Status + Message (optional)             |  "
echo "| ----------------------------------- | ------- | --------------------------------------- |  "

RESULTS_LOG="testResults_$EXENAME/maketest_$date.log"

# copy auxiliary files. All auxiliary files tests/* will be deleted with make clean
cp fchk/*fchk .
cp target/*target .
cp vec/*vec .

#for testfile in `ls H2O*.py`; do
for testfile in `ls *.py`; do

    if [ "$testfile" = "lowdinTestFunctions.py" ]; then
       continue
    fi

    #Run test
    testName=`echo $testfile | gawk '{print substr($1,1,length($1)-3)}'`

    /usr/bin/time -f "%e" -o time.log python3 "$testName.py" "$EXENAME" > output.log 2> error.log 

    status=$( tail -1 output.log )
    error=$( cat error.log ) 
    sed -i '/[a-zA-Z]/d' time.log #remove additional printing...
    duration=$( cat time.log ) 

    printf "| %-35.35s | %-7.7s | %-50.50s |\n" "$testName" "$duration" "$status" | tee -a "$RESULTS_LOG" 

    if [ $(wc -l < output.log) -gt 1 ]; then
    	head -n -1 output.log | tee -a "$RESULTS_LOG" 
    fi

    rm -f output.log 
    rm -f error.log 
    rm -f time.log 

    #Save results 
    find . -maxdepth 1 -name $testName.out -exec mv -t testResults_$EXENAME {} \;
    find . -maxdepth 1 -name $testName"*molden" -exec mv -t testResults_$EXENAME {} \;
    find . -maxdepth 1 -name $testName"*cub" -exec mv -t testResults_$EXENAME {} \;
    find . -maxdepth 1 -name $testName"*dens" -exec mv -t testResults_$EXENAME {} \;
    find . -maxdepth 1 -name $testName"*orb*" -exec mv -t testResults_$EXENAME {} \;
done

failed=$(grep -ch "NOT OK" $RESULTS_LOG)
crashed=$(grep -ch "CRASHED" $RESULTS_LOG)

if [ "$failed" -gt 0 ] ; then     
    echo $failed "tests failed, check output"
    exit $failed
fi
if [ "$crashed" -gt 0 ] ; then     
    echo "$crashed" "tests crashed, check output"
    exit $crashed
fi

if [ "$failed" = 0 ] && [ "$crashed" = 0 ] ; then     
    echo "All tests completed successfully"
    exit 0
fi

