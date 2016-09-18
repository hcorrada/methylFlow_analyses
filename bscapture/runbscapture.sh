#!/bin/bash

### run :  sh sam.sh par1 par2 par3 par4 par5
## Example: sh runbscapture.sh  /Users/faezeh/Projects/methylFlow/data/bscapture ./cpg_output ./region_output 0

### par1 = input file directory

### par2 = cpg output folder name
### par3 = region output folder name




### par4 = 0 > not a sam input
### par4 = 1 > sam input file




## input file is Sam or
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

##input = /cbcb/project-scratch/fdorri/Code/methylFlow/testing/test.sam

####### run with auto lambda ###############################################################

########  evaluation for SAM input ########

pwd=$(pwd)
echo $pwd


datadir=$1
outdir_cpg=$2
outdir_region=$3

conditions="N T"
subjects=$(seq 4 6)
chrs=$(seq 6 9)

for subject in ${subjects}; do
    echo subject: $subject
    for condition in $conditions; do
        echo condition : $conditions
        for chr in $chrs; do
            echo chr : $chrs
            echo "$chr"
            echo "$condition"
            echo "$subject"

            curdir="CAP_${condition}_${subject}"
            echo curdir: $curdir

            if [ ! -d $pwd/$outdir_cpg/$curdir/$chr ];then
                echo "make dir $outdir_cpg/$curdir/${chr}"
                mkdir $pwd/$outdir_cpg/$curdir/${chr}

            fi

            if [ ! -d $pwd/$outdir_region/$curdir/$chr ];then
                echo "make dir $outdir_region/$curdir/${chr}"
                mkdir $pwd/$outdir_region/$curdir/${chr}
            fi


            filename="$datadir/$curdir/${chr}.methylation.withsub.tsv"
            if [ ! -f $filename ]; then
                echo "file $filename does not exist"
                continue
            fi

            echo "run for $filename"



            echo $1
            echo $2
            echo $3
            echo $4
            echo $5



            #export PATH=/Users/faezeh/Libraries/samtools-1.1:$PATH
            #export PATH=/Users/faezeh/Libraries/bcftools-1.1:$PATH

            mf_region="/Users/faezeh/Projects/methylFlow_old/methylFlow/build/methylFlow/methylFlow"

            mf_cpg="${MF_INSTALL_DIR}/bin/methylFlow"

            mfSimulate="${MF_INSTALL_DIR}/bin/mfSimulate"
            mfEvaluate="${MF_INSTALL_DIR}/bin/mfEvaluation"
            avgEvaluate="${MF_INSTALL_DIR}/bin/avgEvaluation"
            samEvaluate="${MF_INSTALL_DIR}/bin/samEvaluation"

            subdirCpG=$outdir_cpg/$curdir/$chr
            subdirRegion=$outdir_region/$curdir/$chr


            echo "Auto lambda - CpG methylflow"
            dirCpG="${pwd}/${subdirCpG}"
            echo $dirCpG

            dirRegion="${pwd}/${subdirRegion}"
            echo $dirRegion

            cd ${dirCpG}

            #cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/sam/auto/
            echo -n "" > methylPercentageRead.txt
            echo -n "" > methylPercentageSam.txt
            echo -n "" > methylPercentageEstimated.txt

            #change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
            #cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/

            if [ "$4" == 1 ]
            then
                echo "SAM input MethylFlow"
                #samtools view -Shu $2 | samtools sort - -o test.sorted | samtools view - -h -o test.sorted.sam
                ${mf_cpg} -i ${filename} -o ${dirCpG} -sam -cpgloss -s 1 -chr ${chr}


            elif [ "$4" == 0 ]
            then
                echo "mf_cpg run"
                ${mf_cpg} -i ${filename} -o ${dirCpG} -cpgloss -s 1 -chr ${chr}
            else

                echo " your input should be 0, 1 "

            fi





            echo "Auto lambda - Region methylFlow"

            echo dirRegion = $dirRegion



            cd ${dirRegion}

            echo -n "" > methylPercentageRead.txt
            echo -n "" > methylPercentageSam.txt
            echo -n "" > methylPercentageEstimated.txt


            if [ "$4" == 1 ]
            then
                echo "SAM input MethylFlow"
                #samtools view -Shu $2 | samtools sort - -o test.sorted | samtools view - -h -o test.sorted.sam
#${mf_region} -i ${filename} -o ${dirRegion} -sam -s 1 chr ${chr}
                ${mf_cpg} -i ${filename} -o ${dirRegion} -sam -s 1 chr ${chr}


            elif [ "$4" == 0 ]
            then
                echo "mf_region run"
#${mf_region} -i ${filename} -o ${dirRegion} -s 1 -chr ${chr}
                ${mf_cpg} -i ${filename} -o ${dirRegion} -s 1 -chr ${chr}

            else

            echo " your input should be 0, 1 "

            fi

        done
    done
done







####### run with non-Auto lambda ###############################################################

