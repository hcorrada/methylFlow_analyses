#! /usr/bin/sh


### this should be change to reading input file name from standard input not hard coded in this!!!


# on CBCB server
# dir <- /cbcb/project-scratch/fdorri/Data/sra
# sh fastq2sam.sh par
## par = input file name
# or
#qsub -q xlarge -l mem=120G,walltime=72:00:00 run.sh -N bismark


input=$1

echo "Hello world"

#######  download sra:
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR124/SRR1248457/SRR1248457.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/wgbs/SRR1248457.sra &


wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR124/SRR1248446/SRR1248446.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/wgbs/SRR1248446.sra &

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR124/SRR1248477/SRR1248477.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/wgbs/SRR1248477.sra &


echo "download cce-2i-120h"


####### sra to fastq:
#/cbcb/project2-scratch/fdorri/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump -A /cbcb/project2-scratch/fdorri/Data/sra/SRR1769243.sra -O /cbcb/project2-scratch/fdorri/Data/sra/fastq &


inputdir="/Users/faezeh/Projects/methylFlow/data/singleCell/fastq/wgbs"
outputdir="/Users/faezeh/Projects/methylFlow/data/singleCell/sam/wgbs"

#echo "start splitting"
#Splitting the fastq files to smaller one so that we can easily run bismark
#split -l 40000000 ${inputdir} /cbcb/project2-scratch/fdorri/Data/sra/fastq/h_
#echo "splitting finished"
cd $inputdir
for f in SRR*.sra
do
echo "$f"
done


#fastq-dump for paired end sequencing
cd $inputdir
for f in SRR*.sra
do
/Users/faezeh/Projects/tools/sratoolkit.2.5.2-mac64/bin/fastq-dump $f -O ${inputdir} &
done
echo "renaming file finished"




#running Bismark and align every fastq files to selected genome using bismark
cd ${inputdir}

im=(${inputdir}/*.fastq)

for ((i=0;i<=${#im[@]};i++))
do
echo "${im1[i]}"
/Users/faezeh/Projects/tools/bismark_v0.14.5/bismark --bowtie2 --non_directional -o ${outputdir}/chr3 /Users/faezeh/Projects/refseq/mouse/chr3 ${im[i]}
done

echo "bismark finished for inputdir1"


#### POOOOOOOOOLING #######


chr=
start=
end=


bam=(${outputdir}/chr3/*.bam)

for ((i=0;i<=${#bam1[@]};i++))
do

samtools sort $bam[i]  ${outputdir}/chr3/sorted.${bam[i]}
echo "samtools sort 1 done"
samtools index ${outputdir}/chr3/sorted.${bam[i]}.bam
samtools view -h -o ${outputdir}/chr3/sorted.${bam[i]}.${start}.bam  ${outputdir}/chr3/sorted.${bam[i]}.bam $chr:$start-$end
echo "number of reads for sorted.${bam1[i]}.${start}.bam"
samtools view -F 0x904 -c ${outputdir}/chr3/sorted.${bam[i]}.${start}.bam

done


#sorting sam files, first need to change to bam and then sort
cd /cbcb/project2-scratch/fdorri/Data/sra/sam/
for f in h_*.sam
do
samtools view -Shu $f | samtools sort -  sorted.${f}
done

echo "sam to sorted bam finished"


#merge all the bam files into a single bam
samtools merge  finall_h.bam sorted.h_*.bam

echo "merging sorted bams finished"


# get final sorted sam from sorted bam
samtools view finall_s.bam -h -o ${outputfile}
echo "sorted bam to final sorted sam file finished"



#### ruunig for another data set #####
#######################################################################################

#Splitting the fastq files to smaller one so that we can easily run bismark
#split -l 40000000 /cbcb/project-scratch/fdorri/Data/sra/fastq/SRR1020523.fastq /cbcb/project-scratch/fdorri/Data/sra/fastq/a_

#renaming the file, adding .fastq extention
#cd /cbcb/project-scratch/fdorri/Data/sra/fastq/
#for f in a_*
#do
#mv $f  ${f}.fastq
#done

#running Bismark and align every fastq files to selected genome using bismark
#cd /cbcb/project-scratch/fdorri/Data/sra/fastq/
#for f in a_*
#do
#/cbcb/project-scratch/fdorri/bismark_v0.11.1/bismark -n 1 -l 50 -o /cbcb/project-scratch/fdorri/Data/sra/sam/ /cbcb/project-scratch/fdorri/Data/refseq/ /cbcb/project-scratch/fdorri/Data/sra/fastq/${f}
#done


#sorting sam files, first need to change to bam and then sort
#cd /cbcb/project-scratch/fdorri/Data/sra/sam/
#for f in b_*.sam
#do
#samtools view -Shu $f | samtools sort -  sorted.${f}
#done


#merge all the bam files into a single bam
#samtools merge  finall_a.bam sorted.a_*.bam

# get final sorted sam from sorted bam
#samtools view finall_a.bam -h -o SRR1020523.sorted.sam

