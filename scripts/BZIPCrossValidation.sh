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

../bin/parana2 pars -t $DATADIR/bZIPData/input/bZIP_n0.adj -u -d $DATADIR/bZIPData/input/bZIP.xml -c -k 20 -p 0.0 -o $OUTDIR/;
