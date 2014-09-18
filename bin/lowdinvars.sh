#!/bin/sh

LOWDIN_HOME="PREFIX/.lowdin2"
LOWDIN_SCRATCH="/scratch"

if [ -z "$LOWDIN_HOME" ]
then
    export LOWDIN_HOME
fi

if [ -z "$LOWDIN_SCRATCH" ]
then
    export LOWDIN_SCRATCH
fi

if [ -z "$LOWDIN_DATA" ]
then
    LOWDIN_DATA="$LOWDIN_HOME/lib"
    export LOWDIN_DATA
fi

if [ -z "$PATH" ]
then
    PATH="$LOWDIN_HOME/bin:$LOWDIN_HOME/utils:"
    export PATH
else
    PATH="${PATH}:$LOWDIN_HOME/bin:$LOWDIN_HOME/utils"
    export PATH
fi

if [ -z "$LIBRARY_PATH" ]
then
    LIBRARY_PATH="$LOWDIN_HOME/lib:$LOWDIN_HOME/utils:"
    export LIBRARY_PATH
else
    LIBRARY_PATH="${LIBRARY_PATH}:$LOWDIN_HOME/bin:$LOWDIN_HOME/utils"
    export LIBRARY_PATH
fi
