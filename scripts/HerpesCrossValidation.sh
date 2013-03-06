#!/bin/bash
# Argument = -d [DATADIR] -o [OUTPUT]

usage()
{
cat << EOF
usage: $0 options

Runs herpes cross validation

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
python ComputeAllEdges.py $DATADIR/HerpesPPIs/OrthologyGroups/OrthologyGroups.proc $DATADIR/HerpesPPIs/trees/rearranged $OUTDIR 60.0

