#! /usr/bin/sh


### this should be change to reading input file name from standard input not hard coded in this!!!


# on CBCB server
# dir <- /cbcb/project-scratch/fdorri/Data/sra
# sh fastq2sam.sh
# or
#qsub -q xlarge -l mem=120G,walltime=72:00:00 run.sh -N bismark


echo "Hello world"

#######  download sra:
#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR101/SRR1015705/SRR1015705.sra -O /cbcb/project-scratch/fdorri/Data/sra/SRR1015705.sra &

####### sra to fastq:
#/cbcb/project-scratch/fdorri/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump -A /cbcb/project-scratch/fdorri/Data/sra/SRR1020509.sra -O /cbcb/project-scratch/fdorri/Data/sra/fastq &



#Splitting the fastq files to smaller one so that we can easily run bismark
 split -l 40000000 /cbcb/project-scratch/fdorri/Data/sra/fastq/SRR1015705.fastq /cbcb/project-scratch/fdorri/Data/sra/fastq/s_



#renaming the file, adding .fastq extention
cd /cbcb/project-scratch/fdorri/Data/sra/fastq/
for f in s_*
do
mv $f  ${f}.fastq
done




#running Bismark and align every fastq files to selected genome using bismark
cd /cbcb/project-scratch/fdorri/Data/sra/fastq/
for f in s_*
do
/cbcb/project-scratch/fdorri/bismark_v0.11.1/bismark -n 1 -l 50 -o /cbcb/project-scratch/fdorri/Data/sra/sam/ /cbcb/project-scratch/fdorri/Data/refseq/ /cbcb/project-scratch/fdorri/Data/sra/fastq/${f}
done

#sorting sam files, first need to change to bam and then sort
cd /cbcb/project-scratch/fdorri/Data/sra/sam/
for f in s_*.sam
do
samtools view -Shu $f | samtools sort -  sorted.${f}
done


#merge all the bam files into a single bam
samtools merge  finall_s.bam sorted.s_*.bam

# get final sorted sam from sorted bam
samtools view finall_s.bam -h -o SRR1015705.sorted.sam



#### ruunig for another data set #####
#######################################################################################

#Splitting the fastq files to smaller one so that we can easily run bismark
split -l 40000000 /cbcb/project-scratch/fdorri/Data/sra/fastq/SRR1020523.fastq /cbcb/project-scratch/fdorri/Data/sra/fastq/a_

#renaming the file, adding .fastq extention
cd /cbcb/project-scratch/fdorri/Data/sra/fastq/
for f in a_*
do
mv $f  ${f}.fastq
done

#running Bismark and align every fastq files to selected genome using bismark
cd /cbcb/project-scratch/fdorri/Data/sra/fastq/
for f in a_*
do
/cbcb/project-scratch/fdorri/bismark_v0.11.1/bismark -n 1 -l 50 -o /cbcb/project-scratch/fdorri/Data/sra/sam/ /cbcb/project-scratch/fdorri/Data/refseq/ /cbcb/project-scratch/fdorri/Data/sra/fastq/${f}
done


#sorting sam files, first need to change to bam and then sort
cd /cbcb/project-scratch/fdorri/Data/sra/sam/
for f in a_*.sam
do
samtools view -Shu $f | samtools sort -  sorted.${f}
done


#merge all the bam files into a single bam
samtools merge  finall_a.bam sorted.a_*.bam

# get final sorted sam from sorted bam
samtools view finall_a.bam -h -o SRR1020523.sorted.sam

