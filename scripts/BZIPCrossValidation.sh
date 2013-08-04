#!/bin/bash
# Argument = -d [DATADIR] -o [OUTPUT]

usage()
{
cat << EOF
usage: $0 options

Runs bZIP cross validation

OPTIONS:
   -h      Show this message
   -d      Data Directory
   -c      Cross validation file
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
         c)  CVFILE=$OPTARG
             ;;
         o)
             OUTDIR=$OPTARG
             ;;
     esac
done

mkdir -p $OUTDIR

../bin/parana2 pars -t $DATADIR/bZIPData/input/bZIP_n0.adj -u -d $DATADIR/bZIPData/input/bZIP.xml -c $CVFILE -k 40 -b 60.0 -p 1.2 -o $OUTDIR/
