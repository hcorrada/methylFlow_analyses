#!/bin/bash

### run :  sh coverage.sh par1 par2

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard



## input file
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> cpgNum >> corrDist ;

# dataFlag = 0  >>> read CpG sites from file(name of the file is at the end of file ( third line))
# dataFlag > 0 >>>> dataFlag equals the number of cpg sites
# dataFlag < 0 >>> read the data from rest of the file

#cpgNum  >>>> Number of cpg or number of lines to be read from file

# freqFlag = 0 >>> randomly choose the frequency of each pattern
# freqFlag = 1 >>> read the frequency of patterns from rest of the file(second line)

pwd=$(pwd)
echo $pwd

mf="${MF_INSTALL_DIR}/bin/methylFlow"
mfSimulate="${MF_INSTALL_DIR}/bin/mfSimulate"
mfEvaluate="${MF_INSTALL_DIR}/bin/mfEvaluation"
avgEvaluate="${MF_INSTALL_DIR}/bin/avgEvaluation"

start=11006910
end=15008000
length=$(($end - $start - 1))



if [ "$1" == 0 ];
then

if [ "$2" == 2 ];
then


echo "Hard Setting"
dir="${pwd}/hard"


cd ${dir}
#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/
########  evaluation for different coverages ########
printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


for i in $(seq 5 3 50)
do
#printf "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/eval.txt
#echo threshold  abdncError  methylCallError TP  FN  FP >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/eval.txt
#printf "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/input.txt

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt
echo 1 $start $length 100 10 1 $i 0 0 100 20 >> ${dir}/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> ${dir}/input.txt
echo "/Users/faezeh/Projects/methylFlow/exps/sam/SRR1015434-11006910-15008000/cpgs.tsv" >> ${dir}/input.txt


#echo 1 757121 1000 70 10 1 $i 0 80 1 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/input.txt
#echo 10 10 10 10 10 10 10 10 10 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/input.txt
#echo $i >> evalCoverage.txt
#printf "   " >> evalCoverage.txt

echo $i
for j in {1..100}
do


#printf "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/shortRead.txt
#printf "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/simPattern.txt
#printf "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/patterns.tsv
#printf "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/weight.txt
#printf "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/match.txt


printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt


#change directory to  /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCoverage"
#../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCoverage"
#../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/  -l 1 -s 1 -chr 1
${mf} -i ${dir}/shortRead.txt -o ${dir} -l 1 -s 1 -chr 1

echo "EvaluateCoverage"
#../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard 757121 757353 $i
${mfEvaluate} ${dir} ${dir} $start $end $i

done
echo "avgEval Start"
#../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard $i
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done



elif [ "$2" == 1 ];
then

echo "Moderate Setting"

dir="${pwd}/moderate"


cd ${dir}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


for i in $(seq 5 3 50)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt
echo 1 $start $length 100 4 1 $i 0 0 100 20 >> ${dir}/input.txt
echo 15 15 35 35 >> ${dir}/input.txt
echo "/Users/faezeh/Projects/methylFlow/exps/sam/SRR1015434-11006910-15008000/cpgs.tsv" >> ${dir}/input.txt


echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt



echo "SimulateCoverage"
${mfSimulate} ${dir}/input.txt ${dir}


echo "MethylFlowCoverage"
${mf} -i ${dir}/shortRead.txt -o ${dir} -s 1 -chr 1

echo "EvaluateCoverage"
${mfEvaluate} ${dir} ${dir} $start $end $i

done
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done




elif [ "$2" == 0 ];
then


echo "Simple Setting"
dir="${pwd}/simple"


cd ${dir}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


for i in $(seq 5 3 50)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt
echo 1 $start $length 100 2 1 $i 0 0 100 20 >> ${dir}/input.txt
echo 25 75 >> ${dir}/input.txt
echo "/Users/faezeh/Projects/methylFlow/exps/sam/SRR1015434-11006910-15008000/cpgs.tsv" >> ${dir}/input.txt


echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt



echo "SimulateCoverage"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCoverage"
${mf} -i ${dir}/shortRead.txt -o ${dir} -s 1 -chr 1

echo "EvaluateCoverage"
${mfEvaluate} ${dir} ${dir} $start $end $i

done
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done



else
echo " your input should be 0 , 1 or 2"
fi



elif [ "$1" == 1 ];
then

if [ "$2" == 2 ];
then

echo "Hard Setting"
dir="${pwd}/hard-Auto"


cd ${dir}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


for i in $(seq 5 3 120)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt
echo 1 757121 1000 70 10 1 $i 0 50 10 >> ${dir}/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> ${dir}/input.txt


echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt



echo "SimulateCoverage"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCoverage"
${mf} -i ${dir}/shortRead.txt -o ${dir} -s 1 -chr 1

echo "EvaluateCoverage"
${mfEvaluate} ${dir} ${dir} 757121 758121 $i

done
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done


elif [ "$2" == 1 ];
then

echo "Moderate Setting"
dir="${pwd}/moderate-Auto"


cd ${dir}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


for i in $(seq 5 3 120)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt
echo 1 757121 1000 70 4 1 $i 0 50 10 >> ${dir}/input.txt
echo 15 15 35 35 >> ${dir}/input.txt


echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt



echo "SimulateCoverage"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCoverage"
${mf} -i ${dir}/shortRead.txt -o ${dir} -s 1 -chr 1

echo "EvaluateCoverage"
${mfEvaluate} ${dir} ${dir} 757121 758121 $i

done
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done


elif [ "$2" == 0 ];
then


echo "Simple Setting"

dir="${pwd}/simple-Auto"


cd ${dir}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


for i in $(seq 5 3 120)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt
echo 1 757121 1000 70 2 1 $i 0 50 10 >> ${dir}/input.txt
echo 25 75 >> ${dir}/input.txt


echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt



echo "SimulateCoverage"
${mfSimulate} ${dir}/input.txt ${dir}


echo "MethylFlowCoverage"
echo ${dir}/shortRead.txt
echo ${dir}
${mf} -i ${dir}/shortRead.txt -o ${dir} -s 1 -chr 1

echo "EvaluateCoverage"
${mfEvaluate} ${dir} ${dir} 757121 758121 $i

done
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done


else
echo " your input should be 0 , 1 or 2"
fi


else
echo " your input should be 0 if it is automatic and 1 if lambda is hard coded"
fi