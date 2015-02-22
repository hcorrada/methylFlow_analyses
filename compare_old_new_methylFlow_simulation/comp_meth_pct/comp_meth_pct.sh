#!/usr/bin/bash

### run :  sh readLength.sh par1

### par1 = 0 > simple
### par1 = 1 > moderate
### par1 = 2 > Hard


## input file
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

# dataFlag = 0  >>> read CpG sites from file
# dataFlag > 0 >>>> dataFlag equals the number of cpg sites
# dataFlag < 0 >>> read the data from rest of the file

# freqFlag = 0 >>> randomly choose the frequency of each pattern
# freqFlag = 1 >>> read the frequency of patterns from rest of the file(second line)

pwd=$(pwd)
echo $pwd

mf_old="/Users/faezeh/Projects/methylFlow_old/methylFlow/build/methylFlow/methylFlow"

mf_new="${MF_INSTALL_DIR}/bin/methylFlow"
mfSimulate="${MF_INSTALL_DIR}/bin/mfSimulate"
mfEvaluate="${MF_INSTALL_DIR}/bin/mfEvaluation"
avgEvaluate="${MF_INSTALL_DIR}/bin/avgEvaluation"






if [ "$1" == 0 ];
then

echo "compare methyl percentage for Simple "
dir_new="${pwd}/Simple/cpg_loss"
dir_old="${pwd}/Simple/region_loss"

echo $dir_new
echo $dir_old

cd ${dir_new}

########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


cd ${dir_old}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt

######

printf "" > ${dir_new}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir_new}/eval.txt

printf "" > ${dir_old}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir_old}/eval.txt


printf "" > ../input.txt

echo 1 757121 3000 70 3 1 20 0 700 10 >> ../input.txt
echo 10 30 60 >> ../input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt


printf "" > ${dir_old}/shortRead.txt
printf "" > ${dir_old}/simPattern.txt
printf "" > ${dir_old}/patterns.tsv
printf "" > ${dir_old}/weight.txt
printf "" > ${dir_old}/match.txt

printf "" > ${dir_new}/shortRead.txt
printf "" > ${dir_new}/simPattern.txt
printf "" > ${dir_new}/patterns.tsv
printf "" > ${dir_new}/weight.txt
printf "" > ${dir_new}/match.txt




echo "SimulateCoverage"
${mfSimulate} ../input.txt ./..


cp ../shortRead.txt ${dir_new}
cp ../shortRead.txt ${dir_old}

cp ../simPattern.txt ${dir_new}
cp ../simPattern.txt ${dir_old}


echo "MethylFlowCoverage_new"
${mf_new} -i ${dir_new}/shortRead.txt -o ${dir_new} -s 1 -chr 1

echo "MethylFlowCoverage_old"
${mf_old} -i ${dir_old}/shortRead.txt -o ${dir_old} -s 1 -chr 1

echo "EvaluateCoverage_new"
${mfEvaluate} ${dir_new} ${dir_new} 757121 757653 $i

echo "EvaluateCoverage_old"
${mfEvaluate} ${dir_old} ${dir_old} 757121 757653 $i

echo "avgEval Start_new"
${avgEvaluate} ${dir_new} ${dir_new} $i

echo "avgEval Start_old"
${avgEvaluate} ${dir_old} ${dir_old} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"






elif [ "$1" == 1 ];
then



echo "compare methyl percentage for Moderate "
dir_new="${pwd}/Moderate/cpg_loss"
dir_old="${pwd}/Moderate/region_loss"

echo $dir_new
echo $dir_old

cd ${dir_new}

########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


cd ${dir_old}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt

######

printf "" > ${dir_new}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir_new}/eval.txt

printf "" > ${dir_old}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir_old}/eval.txt


printf "" > ../input.txt

echo 1 757121 3000 70 10 1 20 0 700 10 >> ../input.txt
echo 5 5 5 5 5 15 15 15 15 15 >> ../input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt


printf "" > ${dir_old}/shortRead.txt
printf "" > ${dir_old}/simPattern.txt
printf "" > ${dir_old}/patterns.tsv
printf "" > ${dir_old}/weight.txt
printf "" > ${dir_old}/match.txt

printf "" > ${dir_new}/shortRead.txt
printf "" > ${dir_new}/simPattern.txt
printf "" > ${dir_new}/patterns.tsv
printf "" > ${dir_new}/weight.txt
printf "" > ${dir_new}/match.txt




echo "SimulateCoverage"
${mfSimulate} ../input.txt ./..


cp ../shortRead.txt ${dir_new}
cp ../shortRead.txt ${dir_old}

cp ../simPattern.txt ${dir_new}
cp ../simPattern.txt ${dir_old}


echo "MethylFlowCoverage_new"
${mf_new} -i ${dir_new}/shortRead.txt -o ${dir_new} -s 1 -chr 1

echo "MethylFlowCoverage_old"
${mf_old} -i ${dir_old}/shortRead.txt -o ${dir_old} -s 1 -chr 1

echo "EvaluateCoverage_new"
${mfEvaluate} ${dir_new} ${dir_new} 757121 757653 $i

echo "EvaluateCoverage_old"
${mfEvaluate} ${dir_old} ${dir_old} 757121 757653 $i

echo "avgEval Start_new"
${avgEvaluate} ${dir_new} ${dir_new} $i

echo "avgEval Start_old"
${avgEvaluate} ${dir_old} ${dir_old} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"






elif [ "$1" == 2 ];
then



echo "compare methyl percentage for Hard "
dir_new="${pwd}/Hard/cpg_loss"
dir_old="${pwd}/Hard/region_loss"

echo $dir_new
echo $dir_old

cd ${dir_new}

########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt


cd ${dir_old}
########  evaluation for different coverages ########

printf "" > evalAvg.txt
echo var'\t'threshold'\t'abdncError'\t'methylCallError'\t'TP'\t'FN'\t'FP >> evalAvg.txt

printf "" > mcf.txt
echo var'\t'minCostFlow >> mcf.txt

printf "" > weight.txt
printf "" > match.txt
printf "" > matchApp.txt

######

printf "" > ${dir_new}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir_new}/eval.txt

printf "" > ${dir_old}/eval.txt
echo threshold  abdncError  methylCallError TP  FN  FP >> ${dir_old}/eval.txt


printf "" > ../input.txt

echo 1 757121 3000 70 30 1 20 0 700 10 >> ../input.txt
echo 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 3 >> ../input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt


printf "" > ${dir_old}/shortRead.txt
printf "" > ${dir_old}/simPattern.txt
printf "" > ${dir_old}/patterns.tsv
printf "" > ${dir_old}/weight.txt
printf "" > ${dir_old}/match.txt

printf "" > ${dir_new}/shortRead.txt
printf "" > ${dir_new}/simPattern.txt
printf "" > ${dir_new}/patterns.tsv
printf "" > ${dir_new}/weight.txt
printf "" > ${dir_new}/match.txt




echo "SimulateCoverage"
${mfSimulate} ../input.txt ./..


cp ../shortRead.txt ${dir_new}
cp ../shortRead.txt ${dir_old}

cp ../simPattern.txt ${dir_new}
cp ../simPattern.txt ${dir_old}


echo "MethylFlowCoverage_new"
${mf_new} -i ${dir_new}/shortRead.txt -o ${dir_new} -s 1 -chr 1

echo "MethylFlowCoverage_old"
${mf_old} -i ${dir_old}/shortRead.txt -o ${dir_old} -s 1 -chr 1

echo "EvaluateCoverage_new"
${mfEvaluate} ${dir_new} ${dir_new} 757121 757653 $i

echo "EvaluateCoverage_old"
${mfEvaluate} ${dir_old} ${dir_old} 757121 757653 $i

echo "avgEval Start_new"
${avgEvaluate} ${dir_new} ${dir_new} $i

echo "avgEval Start_old"
${avgEvaluate} ${dir_old} ${dir_old} $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"






else
echo "Your input should be 0, 1 or 2"
exit 2
fi



