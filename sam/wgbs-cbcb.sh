#! /usr/bin/sh


### this should be change to reading input file name from standard input not hard coded in this!!!


# on CBCB server
# dir <- /cbcb/project-scratch/fdorri/Data/sra
#qsub -q large -l mem=120G,walltime=72:00:00 wgbs.sh -N bismark



refseqdir="/cbcb/project2-scratch/fdorri/Data/refseq/mouse/chr3"
inputdir="/cbcb/project2-scratch/fdorri/Data/wgbs-singleCellBulk/fastq"
outputdir="/cbcb/project2-scratch/fdorri/Data/wgbs-singleCellBulk/sam"

echo "Hello world"

#######  download sra:
#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR124/SRR1248457/SRR1248457.sra -O ${inputdir}/SRR1248457.sra &


#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR124/SRR1248446/SRR1248446.sra -O ${inputdir}/SRR1248446.sra &

#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR124/SRR1248477/SRR1248477.sra -O ${inputdir}/SRR1248477.sra &
#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR124/SRR1248469/SRR1248469.sra -O ${inputdir}/SRR1248469.sra &


#wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR124/SRR1248497/SRR1248497.sra -O ${inputdir}/SRR1248497.sra &



echo "download wgbs single cell"


####### sra to fastq:
#/cbcb/project2-scratch/fdorri/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump -A /cbcb/project2-scratch/fdorri/Data/sra/SRR1769243.sra -O /cbcb/project2-scratch/fdorri/Data/sra/fastq &


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
cd ${inputdir}

im=(${inputdir}/*.sra)
for ((i=0;i<=${#im[@]};i++))
do
/cbcb/project2-scratch/fdorri/sratoolkit.2.5.2-centos_linux64/bin/fastq-dump ${im[i]} -O ${inputdir} &
done
echo "renaming file finished"




#running Bismark and align every fastq files to selected genome using bismark
cd ${inputdir}

im=(${inputdir}/*.fastq)

for ((i=0;i<=${#im[@]};i++))
do
echo "${im[i]}"
/cbcb/project2-scratch/fdorri/bismark_v0.14.5/bismark --bowtie2 --non_directional -o ${outputdir}/chr3 ${refseqdir} ${im[i]}
done

echo "bismark finished for inputdir"



bam=(${outputdir}/chr3/*.bam)

for ((i=0;i<=${#bam[@]};i++))
do

samtools sort ${bam[i]}  ${bam[i]}.sorted
echo "samtools sort $i done"
samtools index ${bam[i]}.sorted.bam

done


#### POOOOOOOOOLING #######


chr=3
start=300000
end=30000000

bam=(${outputdir}/chr3/*.sorted.bam)

for ((i=0;i<${#bam[@]};i++))
do
echo "$i"
#echo "number of reads for ${bam[i]}"
samtools view -h -o ${bam[i]}.${start}.${end}.bam ${bam[i]} ${chr}:${start}-${end}
#echo "number of reads for ${bam[i]}.${start}.${end}.bam"
#samtools view -F 0x904 -c ${bam[i]}.${start}.${end}.bam
samtools view -h -o ${bam[i]}.${start}.${end}.sam ${bam[i]}.${start}.${end}.bam

done



#rm ${outputdir}/chr3/temp.bam

#bam=(${outputdir}/chr3/*.${start}.${end}.bam)
#cp ${bam[0]} ${outputdir}/chr3/out.bam

#for ((i=1;i<${#bam[@]};i++))
#do
#echo "merge bam $i"
#samtools merge ${outputdir}/chr3/temp.bam ${bam[i]} ${outputdir}/chr3/out.bam
#mv ${outputdir}/chr3/temp.bam ${outputdir}/chr3/out.bam

#done



#echo "creating sam files - sam1"
#samtools view -h -o ${outputdir}/chr3/out.sam ${outputdir}/chr3/out.bam


