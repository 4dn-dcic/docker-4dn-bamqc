#!/bin/bash
# Creates bamqc for annotated bam files

INPUT=$1
CHROMSIZES=$2
OUTDIR=$3

FILE_BASE=$(basename $INPUT)
FILE_NAME=${FILE_BASE%.*}

if [ ! -d "$OUTDIR" ]
then
	mkdir $OUTDIR
fi

python3 /usr/local/bin/get_bamqc.py $INPUT $CHROMSIZES $OUTDIR $FILE_NAME
