#!/bin/bash
### example run in two steps
### intermediate file generation for Rscript

### get the absolute path
SPATH="$( cd "$(dirname "$0")" ; pwd -P )"

### ./run_bc.sh P1-noXY.mutect.vcf P1-noXY.battenberg.txt P1-noXY.cellularity_ploidy.txt P1
if [ $# != 4 ]; then
    echo "Please enter: ./run_bc.sh <mutect vcf file> <battenberg file> <cellularity_ploidy file> <out_prefix>"
	exit
else
    echo "running codes ... "
	$SPATH/parseInputData_smc $1 $2 temp.out
	### Rscript run
	Rscript $SPATH/BayCloneC_smc.R temp.out $3 $4 > run1.log 2>&1
fi

