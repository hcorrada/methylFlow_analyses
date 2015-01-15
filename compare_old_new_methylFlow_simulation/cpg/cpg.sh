#!/bin/bash

### run :  sh cpg.sh par1 par2

### par1 = 0 > New mf
### par1 = 1 > Old mf

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

mf_new="/Users/faezeh/Projects/methylFlow/install/bin/methylFlow"
mf_old="/Users/faezeh/Projects/methylFlow_old/methylFlow/build/methylFlow/methylFlow"


mfSimulate="/Users/faezeh/Projects/methylFlow/install/bin/mfSimulate"
mfEvaluate="/Users/faezeh/Projects/methylFlow/install/bin/mfEvaluation"
avgEvaluate="/Users/faezeh/Projects/methylFlow/install/bin/avgEvaluation"



####### run with auto lambda ###############################################################

########  evaluation for different number of CpG sites ########
if [ "$1" == 0 ]
then

if [ "$2" == 2 ]
then

echo "Hard Setting"

dir="${pwd}/hard/new"
cd ${dir}

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 120)
do
printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 70 10 1 20 0 $i 10  >> ${dir}/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> ${dir}/input.txt
#echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCpG"

#../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCpG"
#../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto/  -s 1 -chr 1
${mf_new} -i ${dir}/shortRead.txt -o ${dir} -l 1 -s 1 -chr 1


echo "EvaluateCpG"
#../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto 757121 757353 $i
${mfEvaluate} ${dir} ${dir} 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"
#../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto $i
${avgEvaluate} ${dir} ${dir} $i


echo "avgCpGEval end"

done


elif [ "$2" == 1 ]
then

echo "Moderate Setting"

dir="${pwd}/moderate/new"
cd ${dir}

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 120)
do
printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 70 10 1 20 0 $i 10  >> ${dir}/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> ${dir}/input.txt
#echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCpG"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCpG"
${mf_new} -i ${dir}/shortRead.txt -o ${dir} -l 1 -s 1 -chr 1


echo "EvaluateCpG"
${mfEvaluate} ${dir} ${dir} 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i


echo "avgCpGEval end"
done


elif [ "$2" == 0 ]
then

echo "Simple Setting"

dir="${pwd}/simple/new"
cd ${dir}

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 120)
do
printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 70 10 1 20 0 $i 10  >> ${dir}/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> ${dir}/input.txt
#echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt



echo "SimulateCpG"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCpG"
${mf_new} -i ${dir}/shortRead.txt -o ${dir} -l 1 -s 1 -chr 1


echo "EvaluateCpG"
${mfEvaluate} ${dir} ${dir} 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i


echo "avgCpGEval end"

done


else

echo " your input should be 0, 1 or 2"

fi



####### run with hard coded lambda ###############################################################


########  evaluation for different number of CpG sites ########
elif [ "$1" == 1 ]
then

if [ "$2" == 2 ]
then

echo "Hard Setting"

dir="${pwd}/hard/old"
cd ${dir}

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 120)
do
printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 70 10 1 20 0 $i 10  >> ${dir}/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> ${dir}/input.txt
echo $i

for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt


echo "SimulateCpG"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCpG"
${mf_old} -i ${dir}/shortRead.txt -o ${dir} -l 1 -s 1 -chr 1


echo "EvaluateCpG"
${mfEvaluate} ${dir} ${dir} 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i


echo "avgCpGEval end"

done


elif [ "$2" == 1 ]
then

echo "Moderate Setting"

dir="${pwd}/moderate/old"
cd ${dir}

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 120)
do
printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 70 4 1 20 0 $i 10  >> ${dir}/input.txt
echo 15 15 35 35 >> ${dir}/input.txt
#echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCpG"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCpG"
${mf_old} -i ${dir}/shortRead.txt -o ${dir} -l 1 -s 1 -chr 1


echo "EvaluateCpG"
${mfEvaluate} ${dir} ${dir} 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i


echo "avgCpGEval end"

done


elif [ "$2" == 0 ]
then

echo "Simple Setting"

dir="${pwd}/simple/old"
cd ${dir}

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 120)
do
printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 757121 230 70 2 1 20 0 $i 10  >> ${dir}/input.txt
echo 25 75 >> ${dir}/input.txt
#echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
for j in {1..100}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt


echo "SimulateCpG"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCpG"
${mf_old} -i ${dir}/shortRead.txt -o ${dir} -l 1 -s 1 -chr 1


echo "EvaluateCpG"
${mfEvaluate} ${dir} ${dir} 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i


echo "avgCpGEval end"

done


else

echo " your input should be 0, 1 or 2"

fi


else

echo " your input should be 0 for auto and 1 for hard coded lambda"

fi




















