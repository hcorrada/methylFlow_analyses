#!/bin/bash

### run :  sh sam.sh par1 par2 par3 par4 par5
## Example: sh sam.sh 0 /Users/faezeh/Projects/methylFlow/data/WGS/SRR1020509.sorted.sam SRR1020509 1 3 3100000 15000000

## Example: sh sam.sh 0 /Users/faezeh/Projects/methylFlow/data/singleCell/wgbs/lane7_RSC9N_BS_merged_L007_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_bismark_bt2.deduplicated.bam.sorted.bam.3.3000000.10000000.sam wgbs-bulk 1 3 3020500 3020937
#sh sam.sh 0 /Users/faezeh/Projects/methylFlow/data/bscapture/CAP_N_5/10.methylation.withsub.tsv test 0 3 3000000 3500000
#sh sam.sh 0 /Users/faezeh/Projects/methylFlow/data/WGS/SRR1015434.sorted.sam test1 1 3 3000000 3500000

#${mf} -i /Users/faezeh/Projects/methylFlow/data/singleCell/wgbs/lane8_RSC9N_BI_merged_L008_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_bismark_bt2.deduplicated.bam.sorted.bam.3.3000000.10000000.sam -o lane8 -sam -s 1 -chr 3 -start 3000000 -end 10000000

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

####### run with auto lambda ###############################################################

########  evaluation for SAM input ########


echo $2
echo $3
echo $4
echo $5
echo $6
echo $7


chr=$5
start=$6
end=$7



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



dir="${pwd}/${subdir}"
mkdir ${dir}


samtools view -hb -F 16 $2 > ${dir}/strandf.bam


samtools sort ${dir}/strandf.bam  ${dir}/sorted.1
echo "samtools sort 1 done"


#need indexed bam file to do select a region"

samtools index ${dir}/sorted.1.bam

echo "samtools index done"


samtools view -h -o ${dir}/sorted_region_${start}.1.bam  ${dir}/sorted.1.bam $chr:$start-$end


echo "samtools extract region done"


echo "number of reads for sorted.1"
samtools view -F 0x904 -c ${dir}/sorted_region_${start}.1.bam


file1="${dir}/sorted_region_${start}.1.bam"
samtools view -h -o ${dir}/out.sam $file1


#cd ${dir}

#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/sam/auto/
echo -n "" > methylPercentageRead.txt
echo -n "" > methylPercentageSam.txt
echo -n "" > methylPercentageEstimated.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
if [ "$1" == 0 ]
then
echo "Auto lambda"

if [ "$4" == 1 ]
then
echo "SAM input MethylFlow"

echo 
#samtools view -Shu $2 | samtools sort - -o test.sorted | samtools view - -h -o test.sorted.sam
${mf} -i ${dir}/out.sam -o ${dir} -sam -s 1 -chr $5 -start $6 -end $7
echo "start= $6"

echo "end= $7"


elif [ "$4" == 0 ]
then
${mf} -i $2 -o ${dir} -s 1 -chr $5 -start $6 -end $7
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
${mf} -i $2 -o ${dir} -sam -s 1 -chr $5 -start $6 -end $7

elif ["$3" == 0 ]
then
${mf} -i $2 -o ${dir}  -s 1 -chr $5 -start $6 -end $7

else

echo " your input should be 0, 1 "

fi


else
echo " your input should be 0, 1 : 0 for Auto-Lambda, 1 for hard coded lambda"
fi


