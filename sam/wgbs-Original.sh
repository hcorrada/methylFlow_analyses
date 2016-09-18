#! /usr/bin/sh


### this should be change to reading input file name from standard input not hard coded in this!!!


# on CBCB server
# dir <- /cbcb/project-scratch/fdorri/Data/sra
# sh fastq2sam.sh par
## par = input file name
# or
#qsub -q xlarge -l mem=120G,walltime=72:00:00 run.sh -N bismark


#input=$1

echo "Hello world"

#inputdir="/Users/faezeh/Projects/methylFlow/data/singleCell/fastq/wgbs"
#outputdir="/Users/faezeh/Projects/methylFlow/data/singleCell/sam/wgbs"

inputdir="/cbcb/project2-scratch/fdorri/Data/wgbs-singleCellOriginal/selectedbam"
outputdir="/cbcb/project2-scratch/fdorri/Data/wgbs-singleCellOriginal/sam"

#### POOOOOOOOOLING #######


chr=$1
start=$2
end=$3


bam=(${inputdir}/*.bam)

for ((i=0;i<=${#bam[@]};i++))
do

samtools sort ${bam[i]}  ${inputdir}/chr3/sorted.${bam[i]}
echo "samtools sort 1 done"
samtools index ${inputdir}/chr3/sorted.${bam[i]}.bam
#samtools view -h -o ${outputdir}/chr3/sorted.${bam[i]}.${chr}.${start}.${end}.sam  ${inputdir}/chr3/sorted.${bam[i]}.bam $chr:$start-$end
#echo "number of reads for sorted.${bam1[i]}.${start}.bam"
#samtools view -F 0x904 -c ${outputdir}/chr3/sorted.${bam[i]}.${start}.bam

done


bam=(${inputdir}/sorted.*.bam)

for ((i=0;i<=${#bam[@]};i++))
do

#samtools sort $bam[i]  ${inputdir}/chr3/sorted.${bam[i]}
#echo "samtools sort 1 done"
#samtools index ${inputdir}/chr3/sorted.${bam[i]}.bam
samtools view -h -o ${outputdir}/chr3/${bam[i]}.${chr}.${start}.${end}.sam ${bam[i]} $chr:$start-$end
#echo "number of reads for sorted.${bam1[i]}.${start}.bam"
#samtools view -F 0x904 -c ${outputdir}/chr3/sorted.${bam[i]}.${start}.bam

done


#merge all the bam files into a single bam
#samtools merge  finall_h.bam sorted.h_*.bam

#echo "merging sorted bams finished"


# get final sorted sam from sorted bam
#samtools view finall_s.bam -h -o ${outputfile}
#echo "sorted bam to final sorted sam file finished"

