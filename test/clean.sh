#!/bin/bash

find . -maxdepth 1 -name "*.out" -exec rm -f {} \;
find . -maxdepth 1 -name "*.molden" -exec rm -f {} \;
find . -maxdepth 1 -name "*.cub" -exec rm -f {} \;
find . -maxdepth 1 -name "*.orb*" -exec rm -f {} \;
find . -name "*.vec" -exec rm -f {} \;
find . -name "*.aux" -exec rm -f {} \;
find . -name "*.bas" -exec rm -f {} \;
find . -name "*.sys" -exec rm -f {} \;
find . -name "*.wfn" -exec rm -f {} \;
find . -name "*.pyc" -exec rm -f {} \;
find . -name "*.47" -exec rm -f {} \;
find . -name "*.wfn" -exec rm -f {} \;
find . -name "*.wfx" -exec rm -f {} \;
find . -name "*.eps" -exec rm -f {} \;
find . -name "*.gnp" -exec rm -f {} \;
find . -name "*.vec" -exec rm -f {} \;
find . -name "*.ci" -exec rm -f {} \;
find . -name "*.pyc" -exec rm -f {} \;
find . -name "*.ints" -exec rm -f {} \;
find . -name "*.dens" -exec rm -f {} \;
find . -name "*.coords" -exec rm -f {} \;
find . -name "*.states" -exec rm -f {} \;
