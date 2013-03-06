#!/bin/bash
# Argument = -d [DATADIR] -o [OUTPUT]

usage()
{
cat << EOF
usage: $0 options

Reconstructs the bZIP network at noise levels 0, 10 and 20

OPTIONS:
   -h      Show this message
   -d      Data Directory
   -o      Output Directory
EOF
}

DATATDIR=
OUTDIR=
while getopts “hd:o:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         d)
             DATADIR=$OPTARG
             ;;
         o)
             OUTDIR=$OPTARG
             ;;
     esac
done

mkdir -p $OUTDIR

noiseLevels=( n0 ) #n10 n20 )
penalties=( 1.0 ) #6.0 7.0 )
count=${#noiseLevels[@]}

for i in `seq 1 $count`
do
    n=${noiseLevels[$i-1]}
    p=${penalties[$i-1]}
    echo "../bin/parana2 pars -t ${DATADIR}/bZIPData/input/bZIP_${n}.adj -u -d ${DATADIR}/bZIPData/input/bZIP.xml -k 100 -p ${p} -o $OUTDIR/pars.${n}.edges";
	../bin/parana2 pars -t ${DATADIR}/bZIPData/input/bZIP_${n}.adj -u -d ${DATADIR}/bZIPData/input/bZIP.xml -k 100 -p ${p} -o $OUTDIR/pars.${n}.edges;
done
