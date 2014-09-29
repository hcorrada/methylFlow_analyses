#!/bin/bash

### run :  sh colon.sh par1 par2 par3
## Example: sh colon/colon.sh 0 CAP_N_4 2122  244444444444

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded

### par2 = input folder name

### par2 = start
### par3 = end


## input file is Sam or
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

##input = /cbcb/project-scratch/fdorri/Code/methylFlow/testing/test.sam

####### run with auto lambda ###############################################################

########  evaluation for SAM input ########


echo $2
echo $3
echo $4

if [ "$1" == 0 ]
then
echo "Auto lambda"

for i in `seq 1 22`;
do
echo $i
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/colon/auto/$2/$i/
echo -n "" > methylPercentageRead.txt
echo -n "" > methylPercentageSam.txt
echo -n "" > methylPercentageEstimated.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/ß
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
echo "/cbcb/project-scratch/lmendelo/cancermethylation/methylation/$2/$i.methylation.withsub.tsv"

../build/methylFlow/methylFlow -i /cbcb/project-scratch/lmendelo/cancermethylation/methylation/$2/$i.methylation.withsub.tsv -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/colon/auto/$2/$i/ -s 1 -chr $i


done




fi

#echo "EvaluateCpG"
#../build/samEvaluation/samEvaluation $2 /cbcb/project-scratch/fdorri/Code/methylFlow/testing/sam/auto/patterns.tsv /cbcb/project-scratch/fdorri/Code/methylFlow/testing/sam/auto $4 $5 $3




####### run with non-Auto lambda ###############################################################

if [ "$1" == 1 ]
then
echo "Auto lambda"

for i in `seq 1 22`;
do
echo $i
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/colon/non-Auto/$2/$i/
echo -n "" > methylPercentageRead.txt
echo -n "" > methylPercentageSam.txt
echo -n "" > methylPercentageEstimated.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/ß
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
echo "/cbcb/project-scratch/lmendelo/cancermethylation/methylation/$2/$i.methylation.withsub.tsv"

../build/methylFlow/methylFlow -i /cbcb/project-scratch/lmendelo/cancermethylation/methylation/$2/$i.methylation.withsub.tsv -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/colon/non-Auto/$2/$i/ -s 1 -l 0.5 -chr $i


done


fi

