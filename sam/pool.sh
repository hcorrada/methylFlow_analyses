#! /usr/bin/sh


### this should be change to reading input file name from standard input not hard coded in this!!!


# on CBCB server
# dir <- /cbcb/project-scratch/fdorri/Data/sra
# sh fastq2sam.sh par
## par = input file name
# or
#qsub -q xlarge -l mem=120G,walltime=72:00:00 run.sh -N bismark



#input=$1


mf="${MF_INSTALL_DIR}/bin/methylFlow"

mfSimulate="${MF_INSTALL_DIR}/bin/mfSimulate"
mfEvaluate="${MF_INSTALL_DIR}/bin/mfEvaluation"
avgEvaluate="${MF_INSTALL_DIR}/bin/avgEvaluation"
samEvaluate="${MF_INSTALL_DIR}/bin/samEvaluation"

: <<'end_long_comment'


echo "Hello world"

#######  download sra:
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769139/SRR1769139.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce.4/SRR1769139.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769138/SRR1769138.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce.4/SRR1769138.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769105/SRR1769105.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce.4/SRR1769105.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769104/SRR1769104.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce.4/SRR1769104.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769103/SRR1769103.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce.4/SRR1769103.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769102/SRR1769102.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce.4/SRR1769102.sra &

echo "download cce"


wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769067/SRR1769067.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce_2i_120h.4/SRR1769067.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769069/SRR1769069.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce_2i_120h.4/SRR1769069.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769090/SRR1769090.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce_2i_120h.4/SRR1769090.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769092/SRR1769092.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce_2i_120h.4/SRR1769092.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769094/SRR1769094.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce_2i_120h.4/SRR1769094.sra &
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR176/SRR1769096/SRR1769096.sra -O /Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce_2i_120h.4/SRR1769096.sra &

echo "download cce-2i-120h"
end_long_comment


####### sra to fastq:
#/cbcb/project2-scratch/fdorri/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump -A /cbcb/project2-scratch/fdorri/Data/sra/SRR1769243.sra -O /cbcb/project2-scratch/fdorri/Data/sra/fastq &

outputMethylFlow="/Users/faezeh/Projects/methylFlow/exps/sam"

inputdir1="/Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce_2i_120h.4"
inputdir2="/Users/faezeh/Projects/methylFlow/data/singleCell/fastq/cce.4"

outputdir1="/Users/faezeh/Projects/methylFlow/data/singleCell/sam-single2/cce_2i_120h.4"
outputdir2="/Users/faezeh/Projects/methylFlow/data/singleCell/sam-single2/cce.4"

#echo "start splitting"
#Splitting the fastq files to smaller one so that we can easily run bismark
#split -l 40000000 ${inputdir} /cbcb/project2-scratch/fdorri/Data/sra/fastq/h_
#echo "splitting finished"

: <<'end_long_comment'


cd $inputdir1
for f in SRR*.sra
do
echo "$f"
done


#fastq-dump for paired end sequencing
cd $inputdir1
for f in SRR*.sra
do
/Users/faezeh/Projects/tools/sratoolkit.2.5.2-mac64/bin/fastq-dump -I --split-files $f -O ${inputdir1} &
done
echo "renaming file finished"

cd $inputdir2
for h in SRR*.sra
do
/Users/faezeh/Projects/tools/sratoolkit.2.5.2-mac64/bin/fastq-dump -I --split-files $h -O ${inputdir2} &
done
echo "renaming file finished"


end_long_comment


#running Bismark and align every fastq files to selected genome using bismark
cd ${inputdir1}

im1=(${inputdir1}/*_1.fastq)
im2=(${inputdir1}/*_2.fastq)


: <<'end_long_comment'

for ((i=0;i<=${#im1[@]};i++))
do
echo "${im1[i]}  ,  ${im2[i]}"
/Users/faezeh/Projects/tools/bismark_v0.14.5/bismark  -1 ${im1[i]} -2 ${im2[i]} --bowtie2  -o ${outputdir1}/chr3 /Users/faezeh/Projects/refseq/mouse/chr3
done

end_long_comment



for ((i=0;i<=${#im2[@]};i++))
do
echo "${im1[i]}  ,  ${im2[i]}"
/Users/faezeh/Projects/tools/bismark_v0.14.5/bismark  --bowtie2  -o ${outputdir1}/chr3 /Users/faezeh/Projects/refseq/mouse/chr3 ${im2[i]}
done

echo "bismark finished for inputdir1"




cd ${inputdir2}

mm1=(${inputdir2}/*_1.fastq)
mm2=(${inputdir2}/*_2.fastq)

: <<'end_long_comment'

for ((i=0;i<=${#mm1[@]};i++))
do
echo "${mm1[i]}  ,  ${mm2[i]}"
/Users/faezeh/Projects/tools/bismark_v0.14.5/bismark  -1 ${mm1[i]} -2 ${mm2[i]} --bowtie2  -o ${outputdir2}/chr3 /Users/faezeh/Projects/refseq/mouse/chr3
done
end_long_comment


for ((i=0;i<=${#mm2[@]};i++))
do
echo "${mm1[i]}  ,  ${mm2[i]}"
/Users/faezeh/Projects/tools/bismark_v0.14.5/bismark  --bowtie2  -o ${outputdir2}/chr3 /Users/faezeh/Projects/refseq/mouse/chr3 ${mm2[i]}
done


echo "bismark finished for input dir2"


#### POOOOOOOOOLING #######





bam1=(${outputdir1}/chr3/*.bam)

for ((i=0;i<=${#bam1[@]};i++))
do

samtools sort ${bam1[i]}  ${bam1[i]}.sorted
echo "samtools sort $i done"
samtools index ${bam1[i]}.sorted.bam

done

bam2=(${outputdir2}/chr3/*.bam)

for ((i=0;i<=${#bam1[@]};i++))
do

samtools sort ${bam2[i]}  ${bam2[i]}.sorted
echo "samtools sort $i done"
samtools index ${bam2[i]}.sorted.bam

done

##########################

chr=3
start=300000
end=30000000


bam1=(${outputdir1}/chr3/*.sorted.bam)

for ((i=0;i<${#bam1[@]};i++))
do
echo "$i"
#echo "number of reads for ${bam1[i]}"
samtools view -h -o ${bam1[i]}.${start}.${end}.bam ${bam1[i]} ${chr}:${start}-${end}
echo "number of reads for ${bam1[i]}.${start}.${end}.bam"
samtools view -F 0x904 -c ${bam1[i]}.${start}.${end}.bam

done




bam2=(${outputdir2}/chr3/*.sorted.bam)

for ((i=0;i<${#bam2[@]};i++))
do
samtools view -h -o ${bam2[i]}.${start}.${end}.bam  ${bam2[i]} $chr:$start-$end
echo "number of reads for ${bam1[i]}.${start}.${end}.bam"
samtools view -F 0x904 -c ${bam2[i]}.${start}.${end}.bam

done


rm ${outputdir1}/chr3/temp.bam

bam1=(${outputdir1}/chr3/*.${start}.${end}.bam)
cp ${bam1[0]} ${outputdir1}/chr3/out.bam

for ((i=1;i<${#bam1[@]};i++))
do
echo "merge bam1 $i"
samtools merge ${outputdir1}/chr3/temp.bam ${bam1[i]} ${outputdir1}/chr3/out.bam
mv ${outputdir1}/chr3/temp.bam ${outputdir1}/chr3/out.bam

done




rm ${outputdir2}/chr3/temp.bam

bam2=(${outputdir2}/chr3/*.${start}.${end}.bam)
cp ${bam2[0]} ${outputdir2}/chr3/out.bam

for ((i=1;i<${#bam2[@]};i++))
do
echo "merge bam2 $i"
samtools merge ${outputdir2}/chr3/temp.bam ${bam2[i]} ${outputdir2}/chr3/out.bam
mv ${outputdir2}/chr3/temp.bam ${outputdir2}/chr3/out.bam

done

echo "creating sam files - sam1"
samtools view -h -o ${outputdir1}/chr3/out.sam ${outputdir1}/chr3/out.bam

echo "creating sam files - sam2"
samtools view -h -o ${outputdir2}/chr3/out.sam ${outputdir2}/chr3/out.bam



echo "running methyl flow"
${mf} -i ${outputdir1}/chr3/out.sam -o ${outputMethylFlow}/cce_2i_120h.4_2/chr3 -sam -s 1 -chr $chr -start $start -end $end
${mf} -i ${outputdir2}/chr3/out.sam -o ${outputMethylFlow}/cce.4_2/chr3 -sam -s 1 -chr $chr -start $start -end $end

