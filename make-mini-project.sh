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

# Define region of the alignment which will be used as
# the subset for the mini-project.
region_agtad_mm9="chr11:32000000-33700000"
mini_region=$region_agtad_mm9

# mini-project directory
mp=$PWD

# Get the list of models
mapfile -t model < models.cfg

# loop through each model directory
#     create slurm folder 
#     copy the sample config file 
#     read the list of samples from config
#     loop through each sample in the model
#         create a mini-bam file with the region subset
for m in ${model[@]}; do
    mkdir -p ${m}/slurm
    cd  ${mp}/${m}
    cp "../../${m}/sra_ids.cfg" .
    mapfile -t sra_ids < sra_ids.cfg
    for s in ${sra_ids[@]}; do
        samtools view \
            ../../${m}/${s}.bam \
            $mini_region \
            > ${s}.bam
    done
done