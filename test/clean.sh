#!/bin/bash

if [ -d ./testResults_$1 ]
then
    echo "moving outputs to " ./testResults_$1
    mv *.out ./testResults_$1/
    mv *.molden ./testResults_$1/
    mv *.cub ./testResults_$1/
    mv *.dens ./testResults_$1/
    mv *.orb* ./testResults_$1/
fi

find . -name "*.out" -exec rm -f {} \;
find . -name "*.vec" -exec rm -f {} \;
find . -name "*.aux" -exec rm -f {} \;
find . -name "*.bas" -exec rm -f {} \;
find . -name "*.sys" -exec rm -f {} \;
find . -name "*.wfn" -exec rm -f {} \;
find . -name "*.pyc" -exec rm -f {} \;
find . -name "*.molden" -exec rm -f {} \;
find . -name "*.47" -exec rm -f {} \;
find . -name "*.wfn" -exec rm -f {} \;
find . -name "*.wfx" -exec rm -f {} \;
find . -name "*.orb*" -exec rm -f {} \;
find . -name "*.eps" -exec rm -f {} \;
find . -name "*.gnp" -exec rm -f {} \;
find . -name "*.vec" -exec rm -f {} \;
find . -name "*.ci" -exec rm -f {} \;
find . -name "*.pyc" -exec rm -f {} \;
find . -name "*.ints" -exec rm -f {} \;
find . -name "*.dens" -exec rm -f {} \;
find . -name "*.cub" -exec rm -f {} \;
find . -name "*.coords" -exec rm -f {} \;
find . -name "*.states" -exec rm -f {} \;
