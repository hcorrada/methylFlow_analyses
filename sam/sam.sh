#!/bin/bash

### run :  sh sam.sh par1 par2 par3 par4 par5
## Example: sh sam.sh 0 /Users/faezeh/Projects/methylFlow/data/sam/SRR1020537.sorted.sam SRR1020537-3006910-15008000 1 3006910 15008000

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded



### par2 = input file name
### par3 = output folder name



### par4 = 0 > not a sam input
### par4 = 1 > sam input file


### par5 = start
### par6 = end


## input file is Sam or
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

##input = /cbcb/project-scratch/fdorri/Code/methylFlow/testing/test.sam

####### run with auto lambda ###############################################################

########  evaluation for SAM input ########


echo $2
echo $3
echo $4
echo $5
echo $6


pwd=$(pwd)
echo $pwd

#export PATH=/Users/faezeh/Libraries/samtools-1.1:$PATH
#export PATH=/Users/faezeh/Libraries/bcftools-1.1:$PATH

mf="${MF_INSTALL_DIR}/bin/methylFlow"

mfSimulate="${MF_INSTALL_DIR}/bin/mfSimulate"
mfEvaluate="${MF_INSTALL_DIR}/bin/mfEvaluation"
avgEvaluate="${MF_INSTALL_DIR}/bin/avgEvaluation"
samEvaluate="${MF_INSTALL_DIR}/bin/samEvaluation"

subdir=$3


if [ "$1" == 0 ]
then
echo "Auto lambda"
dir="${pwd}/${subdir}"
mkdir ${dir}


cd ${dir}

#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/sam/auto/
echo -n "" > methylPercentageRead.txt
echo -n "" > methylPercentageSam.txt
echo -n "" > methylPercentageEstimated.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/

if [ "$4" == 1 ]
then
echo "SAM input MethylFlow"
#samtools view -Shu $2 | samtools sort - -o test.sorted | samtools view - -h -o test.sorted.sam
${mf} -i $2 -o ${dir} -sam -s 1 -chr 1 -start $5 -end $6


elif [ "$4" == 0 ]
then
${mf} -i $2 -o ${dir} -s 1 -chr 1 -start $5 -end $6
else

echo " your input should be 0, 1 "

fi





####### run with non-Auto lambda ###############################################################

elif [ "$1" == 1 ]
then



echo "Non-Auto lambda"
echo "Auto lambda"

dir="${pwd}/non-auto/${subdir}"
mkdir ${dir}

cd ${dir}

#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/sam/auto/
echo -n "" > methylPercentageRead.txt
echo -n "" > methylPercentageSam.txt
echo -n "" > methylPercentageEstimated.txt



if [ "$4" == 1 ]
then
echo "SAM input MethylFlow"
#samtools view -Shu $2 | samtools sort - -o test.sorted | samtools view - -h -o test.sorted.sam
${mf} -i $2 -o ${dir} -sam -s 1 -chr 1 -start $5 -end $6

elif ["$3" == 0 ]
then
${mf} -i $2 -o ${dir}  -s 1 -chr 1 -start $5 -end $6

else

echo " your input should be 0, 1 "

fi


else
echo " your input should be 0, 1 : 0 for Auto-Lambda, 1 for hard coded lambda"
fi


