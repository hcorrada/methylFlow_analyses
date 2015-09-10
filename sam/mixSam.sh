#!/bin/bash

### run :  sh sam.sh par1 par2 par3 par4 par5
## Example: sh mixSam.sh /cbcb/project2-scratch/fdorri/Data/sra/sam/mouse/chr3/test1.bam /cbcb/project2-scratch/fdorri/Data/sra/sam/mouse/chr3/test2.bam  test 1 4 3 30 150000000

#sh mixSam.sh ~/Projects/methylFlow/data/singleCell/sam/sorted.SRR1769174_1.fastq_bismark_bt2_pe.bam ~/Projects/methylFlow/data/singleCell/sam/sorted.SRR1769202_1.fastq_bismark_bt2_pe.bam bamMix 1 4 3 3000 15000000

#sh mixSam.sh /Users/faezeh/Projects/methylFlow/data/singleCell/sam/cce_2i_120h.4/chr3/out.bam /Users/faezeh/Projects/methylFlow/data/singleCell/sam/cce.4/chr3/out.bam  poolMix 1 1 3 3019574 3019795

#sh mixSam.sh /Users/faezeh/Projects/methylFlow/data/singleCell/sam/wgbs/chr3/SRR1248477.fastq_bismark_bt2.bam /Users/faezeh/Projects/methylFlow/data/singleCell/sam/wgbs/chr3/SRR1248457.fastq_bismark_bt2.bam  wgbsMix57-77 1 1 3 30000 30000000


### par1 = input file name 1
### par2 = input file name 2
### par3 = output folder name



### par4 =  input file 1 ratio
### par5 =  input file 2 ratio

### par6 = chr
### par7 = start
### par8 = end


## input file is Sam or
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

##input = /cbcb/project-scratch/fdorri/Code/methylFlow/testing/test.sam

####### run with auto lambda ###############################################################

########  evaluation for SAM input ########

echo $1
echo $2
echo $3
echo $4
echo $5
echo $6
echo $7
echo $8


pwd=$(pwd)
echo $pwd

#export PATH=/Users/faezeh/Libraries/samtools-1.1:$PATH
#export PATH=/Users/faezeh/Libraries/bcftools-1.1:$PATH

mf="${MF_INSTALL_DIR}/bin/methylFlow"

mfSimulate="${MF_INSTALL_DIR}/bin/mfSimulate"
mfEvaluate="${MF_INSTALL_DIR}/bin/mfEvaluation"
avgEvaluate="${MF_INSTALL_DIR}/bin/avgEvaluation"
samEvaluate="${MF_INSTALL_DIR}/bin/samEvaluation"


ratio1=$4
ratio2=$5

chr=$6
start=$7
end=$8

subdir=$3


dir="${pwd}/${subdir}"
if [ ! -f ${dir} ]
then
mkdir ${dir}
fi
#samtools view $1 chr${chr}:${start}-${end} > in1.bam

#samtools view $2 chr${chr}:${start}-${end} > in2.bam

samtools sort $1  ${dir}/sorted.1
echo "samtools sort 1 done"
samtools sort $2  ${dir}/sorted.2
echo "samtools sort 2 done"


#need indexed bam file to do select a region"

samtools index ${dir}/sorted.1.bam
samtools index ${dir}/sorted.2.bam

echo "samtools index done"


samtools view -h -o ${dir}/sorted_region_${start}.1.bam  ${dir}/sorted.1.bam $chr:$start-$end


samtools view -h -o ${dir}/sorted_region_${start}.2.bam  ${dir}/sorted.2.bam $chr:$start-$end

echo "samtools extract region done"


echo "number of reads for sorted.1"
samtools view -F 0x904 -c ${dir}/sorted_region_${start}.1.bam


echo "number of reads for sorted.2"
samtools view -F 0x904 -c ${dir}/sorted_region_${start}.2.bam


file1="${dir}/sorted_region_${start}.1.bam"
file2="${dir}/sorted_region_${start}.2.bam"

#: <<'end_long_comment'

if [ "$ratio1" == 1 ]
then
for i in $(seq 1 1 $ratio2)
do
if [ ! -f ${dir}/out.bam ]
then
cp $file1 ${dir}/out.bam
fi

echo "$i"
echo "ratio1=1"
mv ${dir}/out.bam ${dir}/temp.bam
samtools merge ${dir}/out.bam $file2 ${dir}/temp.bam


done

fi

if [ "$ratio2" == 1 ]
then
for i in $(seq 1 1 $ratio1)
do
if [ ! -f ${dir}/out.bam ]
then
cp $file2 ${dir}/out.bam
fi

echo "$i"
echo "ratio2=1"
mv ${dir}/out.bam ${dir}/temp.bam
samtools merge ${dir}/out.bam $file1 ${dir}/temp.bam

done

fi


samtools sort ${dir}/out.bam  ${dir}/out.sorted

samtools view -h -o ${dir}/out.sam ${dir}/out.sorted.bam

samtools view -h -o ${dir}/out1.sam $file1
samtools view -h -o ${dir}/out2.sam $file2

if [ ! -f ${dir}/1 ]
then
mkdir ${dir}/1
fi

if [ ! -f ${dir}/2 ]
then
mkdir ${dir}/2
fi

cd ${dir}

#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/sam/auto/
echo -n "" > methylPercentageRead.txt
echo -n "" > methylPercentageSam.txt
echo -n "" > methylPercentageEstimated.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
#cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SAM input MethylFlow"
#samtools view -Shu $2 | samtools sort - -o test.sorted | samtools view - -h -o test.sorted.sam

end=$((end + 100))

${mf} -i ${dir}/out.sam -o ${dir} -sam -s 1 -chr $chr -start $start -end $end
${mf} -i ${dir}/out1.sam -o ${dir}/1 -sam -s 1 -chr $chr -start $start -end $end
${mf} -i ${dir}/out2.sam -o ${dir}/2 -sam -s 1 -chr $chr -start $start -end $end

#end_long_comment
