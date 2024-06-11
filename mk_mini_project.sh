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

region_agtad_mm9="chr11:32000000-33700000"
mini_region=$region_agtad_mm9

mapfile -t model < models.cfg


for i in ${!model[@]}; do
    mkdir -p ${model[$i]}/slurm
done

for i in ${!model[@]}; do
    samtools view \
        ../${model[i]}/${model[i]}.bam \
        $mini_region > \
        ${model[i]}/${model[i]}.bam
done