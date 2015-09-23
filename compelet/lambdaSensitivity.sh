#!/usr/bin/bash

### run :  sh readLength.sh par1 par2



### par1 = 0 > simple
### par1 = 1 > moderate
### par1 = 2 > Hard


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


if [ "$1" == 2 ];
then


echo "Hard Setting"
dir="${pwd}/hard-lambda"
if [ ! -f ${dir} ]
then
mkdir ${dir}
fi

cd ${dir}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


i=0.03125
i=0.001
for k in `seq 0 1 18`
do
i=$(echo $i*2 | bc)


printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 $start $length 100 10 1 20 0 0 100 20 >> ${dir}/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> ${dir}/input.txt
echo "/Users/faezeh/Projects/methylFlow/exps/sam/SRR1015434-11006910-15008000/cpgs.tsv" >> ${dir}/input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt



for j in {1..10}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt



echo "SimulateCoverage"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCoverage"
${mf} -i ${dir}/shortRead.txt -o ${dir} -s 1 -l $i -chr 1

echo "EvaluateCoverage"
${mfEvaluate} ${dir} ${dir} $start $end $i

done
echo "avgEval Start $i"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"

done

elif [ "$1" == 1 ];
then


echo "Moderate Setting"
dir="${pwd}/moderate-lambda"
if [ ! -f ${dir} ]
then
mkdir ${dir}
fi

cd ${dir}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


#i=0.03125
i=0.001
for k in `seq 0 1 18`
do
i=$(echo $i*2 | bc)


printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 $start $length 100 4 1 20 0 0 100 20 >> ${dir}/input.txt
echo 15 15 35 35 >> ${dir}/input.txt
echo "/Users/faezeh/Projects/methylFlow/exps/sam/SRR1015434-11006910-15008000/cpgs.tsv" >> ${dir}/input.txt

#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt

for j in {1..10}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt



echo "SimulateCoverage"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCoverage"
${mf} -i ${dir}/shortRead.txt -o ${dir} -s 1 -l $i -chr 1

echo "EvaluateCoverage"
${mfEvaluate} ${dir} ${dir} $start $end $i

done
echo "avgEval Start"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"

done





elif [ "$1" == 0 ];
then


echo "Simple Setting"

dir="${pwd}/simple-lambda"
if [ ! -f ${dir} ]
then
mkdir ${dir}
fi

cd ${dir}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt




#i=0.03125
i=0.001
for k in `seq 0 1 18`
do
i=$(echo $i*2 | bc)


printf "" > ${dir}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir}/eval.txt
printf "" > ${dir}/input.txt

echo 1 $start $length 100 2 1 20 0 0 100 20 >> ${dir}/input.txt
echo 25 75 >> ${dir}/input.txt
echo "/Users/faezeh/Projects/methylFlow/exps/sam/SRR1015434-11006910-15008000/cpgs.tsv" >> ${dir}/input.txt

#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt

for j in {1..10}
do

printf "" > ${dir}/shortRead.txt
printf "" > ${dir}/simPattern.txt
printf "" > ${dir}/patterns.tsv
printf "" > ${dir}/weight.txt
printf "" > ${dir}/match.txt



echo "SimulateCoverage"
${mfSimulate} ${dir}/input.txt ${dir}

echo "MethylFlowCoverage"
${mf} -i ${dir}/shortRead.txt -o ${dir} -s 1 -l $i -chr 1

echo "EvaluateCoverage"
${mfEvaluate} ${dir} ${dir} $start $end $i

done
echo "avgEval Start $i"
${avgEvaluate} ${dir} ${dir} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"

done



else
echo "Your input should be 0, 1 or 2"
exit 2
fi

