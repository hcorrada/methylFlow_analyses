#! /bin/bash

datadir=~/willow/project-scratch/lmendelo/cancermethylation/methylation
outdir=reads
conditions="N T"
subjects=$(seq 4 6)
chrs=$(seq 1 22)

for subject in ${subjects}; do
    echo subject: $subject
    for condition in $conditions; do
	curdir="CAP_${condition}_${subject}"
	echo curdir: $curdir

	if [ ! -d $outdir/$curdir ]; then
	    echo "make dir $outdir/$curdir"
	    command=( "mkdir" "-p" "$outdir/$curdir" )
	    echo "${command[@]}"
	    "${command[@]}"
	fi

	for chr in $chrs; do
	    filename="$datadir/$curdir/${chr}.methylation.withsub.tsv"
	    if [ ! -f $filename ]; then
		echo "file $filename does not exist"
		continue
	    fi
	    echo "copy $filename"
	    command=( "cp"  "$filename" "$outdir/$curdir" )
	    echo "${command[@]}"
	    "${command[@]}"
	done
    done
done
