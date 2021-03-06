#!/bin/bash

### run :  sh runmf.sh 

# see ../README.md
#
# assumes methylFlow binary is installed
# in ${MF_INSTALL_DIR}/bin"

# the input data is located at
# "./reads"

# and organized as
# | | | CAP_<condition>_<subject> (for each subject and condition [N|T])
# | | | | <chr>.methylation.withsub.tsv (this is the read file for one chromosme)

# output will be written to 
# "./cpg_output"

# and organized as
# | | | CAP_<condition>_<subject> (output directory for each subject and condition)
# | | | | <chr> (this is where the output for a given chromosome will be placed)

mf="${MF_INSTALL_DIR}/bin/methylFlow"
indir="./reads"
outdir="./cpg_output"

subjects=$(seq 4 6)
conditions="N T"
#chrs=$(seq 1 22)
chrs=$(seq 13 13)

# check command
#${mf} -h && { echo 'command ${mf} -h failed'; exit 1; }

echo "command worked, continuing..."

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

for subject in $subjects; do
    echo "subject: $subject"

    for condition in $conditions; do
	echo "subject: $subject, condition: $condition"

	curindir="${indir}/CAP_${condition}_${subject}"
	if [ ! -d $curindir ]; then
	    echo "directory $curindir not found"
	    continue
	fi

	curoutdir="${outdir}/CAP_${condition}_${subject}"
	echo $curoutdir

	if [ ! -d $curoutdir ]; then
	    echo "making directory ${curoutdir}"
	    mkdir -p $curoutdir
	fi

	for chr in $chrs; do
	    curinfile="${curindir}/${chr}.methylation.withsub.tsv"
	    if [ ! -f $curinfile ]; then
		echo "Input file ${curinfile} not found"
		continue
	    fi

	    tmpoutdir="${curoutdir}/${chr}"
	    if [ ! -d ${tmpoutdir} ]; then
		mkdir -p ${tmpoutdir}
	    fi
	    command=( "$mf" "-i" "${curinfile}" "-o" "${tmpoutdir}" "-s" "1" "-chr" "chr${chr}" "--pctselect" "-e 0.01" )
	    echo "${command[@]}"
	    "${command[@]}"
	done
    done
done

