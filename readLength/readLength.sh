#!/usr/bin/bash

### run :  sh readLength.sh par1 par2

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard


## input file
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

# dataFlag = 0  >>> read CpG sites from file
# dataFlag > 0 >>>> dataFlag equals the number of cpg sites
# dataFlag < 0 >>> read the data from rest of the file

# freqFlag = 0 >>> randomly choose the frequency of each pattern
# freqFlag = 1 >>> read the frequency of patterns from rest of the file(second line)

pwd=$(pwd)
echo $pwd

mf="${MF_INSTALL_DIR}/bin/methylFlow"
mfSimulate="${MF_INSTALL_DIR}/bin/mfSimulate"
mfEvaluate="${MF_INSTALL_DIR}/bin/mfEvaluation"
avgEvaluate="${MF_INSTALL_DIR}/bin/avgEvaluation"

if [ "$1" == 0 ];
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


for i in $(seq 5 3 230)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 $i 10 1 20 0 80 10 >> ${dir}/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> ${dir}/input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt




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
${mfEvaluate} ${dir} ${dir} 757121 757653 $i

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


for i in $(seq 5 3 230)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 $i 4 1 20 0 80 10 >> ${dir}/input.txt
echo 15 15 35 35 >> ${dir}/input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt




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
${mfEvaluate} ${dir} ${dir} 757121 757653 $i

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


for i in $(seq 5 3 230)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 $i 2 1 20 0 80 10 >> ${dir}/input.txt
echo 25 75 >> ${dir}/input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt




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
${mfEvaluate} ${dir} ${dir} 757121 757653 $i

done
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done

else
echo "Your input should be 0, 1 or 2"
exit 2
fi

############################################################# hard coded lambda #####################################

elif [ "$1" == 1 ];
then

if [ "$2" == 2 ];
then


echo "Hard Setting"
dir="${pwd}/hard"


cd ${dir}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


for i in $(seq 5 3 230)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 $i 10 1 20 0 80 10 >> ${dir}/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> ${dir}/input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt




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
${mfEvaluate} ${dir} ${dir} 757121 757653 $i

done
echo "avgEval Start"
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


for i in $(seq 5 3 230)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 $i 4 1 20 0 80 10 >> ${dir}/input.txt
echo 15 15 35 35 >> ${dir}/input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt




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
${mfEvaluate} ${dir} ${dir} 757121 757653 $i

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


for i in $(seq 5 3 230)
do

printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 $i 2 1 20 0 80 10 >> ${dir}/input.txt
echo 25 75 >> ${dir}/input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt




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
${mfEvaluate} ${dir} ${dir} 757121 757653 $i

done
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done




else
echo "Your input should be 0, 1 or 2"
exit 2
fi

else

echo "Your input should be 0 for auto 1 for hard coded lambda"

fi







