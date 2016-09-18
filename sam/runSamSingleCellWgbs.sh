#!/bin/bash

inputdir="/Users/faezeh/Projects/methylFlow/data/singleCell/wgbs"
outputdir="/singleCellwgbs_Ser"


sam=(${inputdir}/lane6_*.sam)


for ((i=0;i<${#sam[@]};i++))
do


sh sam.sh 0 ${sam[i]} ${outputdir}/$i 1 3 3000000 10000000
done

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded



### par2 = input file name
### par3 = output folder name



### par4 = 0 > not a sam input
### par4 = 1 > sam input file

### par5 = chr
### par6 = start
### par7 = end


## input file is Sam or
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

##input = /cbcb/project-scratch/fdorri/Code/methylFlow/testing/test.sam
