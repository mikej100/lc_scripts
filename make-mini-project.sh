#!/bin/bash
#
# Make a mini-project in a subdirectory for testing
# 
# Utilisation
# In the main project create a subdirectory for the mini project
# cd into the mini-project 
# create a file called models.cfg and enter the names of the model
## folders in the project directory.
# Run this script.
# Subdirectoris are created in the mini project for each model
# and put a mini-bam file there containing just the region 
# defined by the mini_region

module load samtools/1.17

region_agtad_mm9="chr11:32000000-33700000"
mini_region=$region_agtad_mm9

mapfile -t model < models.cfg

mp=$PWD
echo "location of mini-project $mp"

for i in ${!model[@]}; do
    mkdir -p ${model[$i]}/slurm
done

for m in ${model[@]}; do
    echo "working on mini-model $m"
    cd  ${mp}/${m}
    echo "Moved to directory ${PWD##}"
    cp "../../${m}/sra_ids.cfg" .
    

    mapfile -t sra_ids < sra_ids.cfg
    
    for s in ${sra_ids[@]}; do
        samtools view \
            ../../${m}/${s}.bam \
            $mini_region \
            > ${s}.bam
    done
done