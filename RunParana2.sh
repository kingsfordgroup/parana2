#!/bin/sh

PLATFORM=$(uname)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $PLATFORM == "Darwin" ]
then
    DYLD_LIBRARY_PATH=${DIR}/../lib:$DYLD_LIBRARY_PATH ${DIR}/parana2 $@
elif [ $PLATFORM == "Linux" ]
then
    LD_LIBRARY_PATH=${DIR}/../lib:$LD_LIBRARY_PATH ${DIR}/parana2 $@
fi
