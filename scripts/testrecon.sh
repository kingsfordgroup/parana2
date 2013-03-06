#!/bin/bash
# Argument = -d [DATADIR] -o [OUTPUT]

ofname=`tempfile`
nl=n$1
p=
case nl in
    "n0" )
        p=1.0;;
    "n10" )
        p=6.0;;
    "n20" )
        p=7.0;;
esac

DATADIR=/Users/rob/SoftwareStaging/Parana2Data/

ofile=`mktemp`

echo "../bin/parana2 pars -t ${DATADIR}/bZIPData/input/bZIP_${n}.adj -u -d ${DATADIR}/bZIPData/input/bZIP.xml -k 100 -p ${p} -o ${ofile}";
../bin/parana2 pars -t ${DATADIR}/bZIPData/input/bZIP_${n}.adj -u -d ${DATADIR}/bZIPData/input/bZIP.xml -k 100 -p ${p} -o ${ofile};
