#!/bin/bash

### run :  sh runmf.sh 

# assumes that directory is laid out as:
# methylFlow
# | tools/bin/methylFlow (this is the executable)
# | src/methylFlow (this is the source directory)
# | exps/bscapture (this script is in this directory)
# | | reads (this is where the data is)
# | | | CAP_<condition>_<subject> (for each subject and condition [N|T])
# | | | | <chr>.methylation.withsub.tsv (this is the read file for one chromosme)
# | | output (this is where the output goes, this directory is created by this script if it doesn't exist)
# | | | CAP_<condition>_<subject> (output directory for each subject and condition)
# | | | | <chr> (this is where the output for a given chromosome will be placed)

# to compile the exectuable so it is found where expected do
# pushd ../../src/methylFlow
# mkdir build && pushd build
# cmake -DCMAKE_INSTALL_PREFIX="../../tools" ..
# make
# make install
# popd && popd

mf="../../tools/bin/methylFlow"
indir="reads"
outdir="cpg_output"

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

