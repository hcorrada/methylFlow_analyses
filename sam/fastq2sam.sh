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
#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR101/SRR1015705/SRR1015705.sra -O /cbcb/project2-scratch/fdorri/Data/sra/SRR1015705.sra &

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769060/SRR1769060.sra -O /cbcb/project2-scratch/fdorri/Data/sra/SRR1769060.sra

####### sra to fastq:
#paired-end
#/cbcb/project2-scratch/fdorri/sratoolkit.2.5.2-centos_linux64/bin/fastq-dump -I --split-files /cbcb/project2-scratch/fdorri/Data/sra/SRR1769245.sra -O /cbcb/project2-scratch/fdorri/Data/sra/fastq &
#single
#/cbcb/project2-scratch/fdorri/sratoolkit.2.5.2-centos_linux64/bin/fastq-dump /cbcb/project2-scratch/fdorri/Data/sra/SRR1769060.sra -O /cbcb/project2-scratch/fdorri/Data/sra/fastq &


inputdir="/cbcb/project-scratch/fdorri/Data/sra/fastq/${input}"
outputfile="${input}.sorted.sam"

echo "start splitting"
#Splitting the fastq files to smaller one so that we can easily run bismark
split -l 40000000 ${inputdir} /cbcb/project2-scratch/fdorri/Data/sra/fastq/h_
echo "splitting finished"



#renaming the file, adding .fastq extention
cd /cbcb/project2-scratch/fdorri/Data/sra/fastq/
for f in h_*
do
mv $f  ${f}.fastq
done
echo "renaming file finished"


#making ref for bismark, seperatly for each chromosome for easier sorting of sam/bam file
/cbcb/project2-scratch/fdorri/bismark_v0.14.5/bismark_genome_preparation --path_to_bowtie /fs/szattic-asmg6/rob/dist/bin/  --verbose /cbcb/project2-scratch/fdorri/Data/refseq/mouse/chr3/ &
# /fs/szattic-asmg6/rob/dist/bin/bowtie2 use bowtie 2

#running Bismark and align every fastq files to selected genome using bismark
cd /cbcb/project2-scratch/fdorri/Data/sra/fastq/
for f in h_*
do

#/cbcb/project2-scratch/fdorri/bismark_v0.11.1/bismark -n 1 -l 50 -o /cbcb/project2-scratch/fdorri/Data/sra/sam/ /cbcb/project2-scratch/fdorri/Data/refseq/ /cbcb/project2-scratch/fdorri/Data/sra/fastq/${f}

#paired-end
/cbcb/project2-scratch/fdorri/bismark_v0.14.5/bismark  -1 /cbcb/project2-scratch/fdorri/Data/sra/fastq/SRR1769245_1.fastq -2 /cbcb/project2-scratch/fdorri/Data/sra/fastq/SRR1769245_2.fastq --bowtie2  -o /cbcb/project2-scratch/fdorri/Data/sra/sam/mouse/chr10 /cbcb/project2-scratch/fdorri/Data/refseq/mouse/chr10
#single
/cbcb/project2-scratch/fdorri/bismark_v0.14.5/bismark --bowtie2 -o /cbcb/project2-scratch/fdorri/Data/sra/sam/mouse/chr10 /cbcb/project2-scratch/fdorri/Data/refseq/mouse/chr10 /cbcb/project2-scratch/fdorri/Data/sra/fastq/SRR1769060.fastq

done

echo "bismark finished"


#if output is sam, then for sorting sam files, first need to change to bam and then sort
cd /cbcb/project2-scratch/fdorri/Data/sra/sam/
for f in h_*.sam
do
samtools view -Shu $f | samtools sort -  sorted.${f}
done

#if output is bam, then you only need to sort
cd /cbcb/project2-scratch/fdorri/Data/sra/sam/mouse/chr3
for f in *.bam
do
samtools sort $f  sorted.${f}
done

echo " sorted bam finished"


#multiple bam#####merge all the bam files into a single bam
samtools merge  finall_h.bam sorted.h_*.bam

echo "merge sorted bams finished"


# get final sorted sam from sorted bam
cd /cbcb/project2-scratch/fdorri/Data/sra/sam/mouse/chr3
for f in sorted.*.bam
do
samtools view $f -h -o ${f}.sam
done


###multiple bam, then sam
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

