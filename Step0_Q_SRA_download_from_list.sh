#!/bin/bash
#SBATCH --partition=long
#SBATCH --output=slurm/%j_%x.out 
#SBATCH --error=slurm/%j_%x.out
#SBATCH --ntasks=4
#SBATCH --mem=12G
#SBATCH --time=00-10:00:00

# Download SRA files specified by list in current folder
# ======================================================
#
# Instructions for use
# Create a directory for the set of files to be processed together (the 'model')
# in the model directory create a file named "sra_ids.cfg" 
# and enter the SRA ids one per line.
# from the model directory enter:
#   "sbatch <path-to-sript-folder>/SRA_download_from_list.sh"
# 
# TODO download file description info so user can confirm is what expected
# Mike Jennings, based on script by Lucy Cornell
#


now() {
    date +"%Y-%m-%dT%T"
}

echo $(now) Starting $(basename "${BASH_SOURCE}")

# Read SRA_ids from a a config file. One per line.
mapfile -t sra_ids < sra_ids.cfg

#echo "Show git commit reference:"
#git show -s --format=%h%x09%ci%x09%an%x09%s

#Model name is the folder name
model=$(echo $(pwd) | awk -F/ '{print $NF}') ;

#Create a date and timestamp for when analysis began
Start_time=`date`
Time=`date +"%T"`

module load sratoolkit
# vdb-config --interactive

now() {
    date +"%Y-%m-%dT%T"
}

echo $(now) Starting $(basename "${BASH_SOURCE}")
# Show git info for scripts folder
${SCRIPTS}/scripts_info.sh || true

#TODO Understand requirement for creating meaningful filename
echo "SRA ids: "${sra_ids[@]}
for sra_id in ${sra_ids[@]}; do
    fastq-dump --split-files $sra_id 
    mv $sra_id"_1.fastq"  $sra_id"_R1.fastq"
    mv $sra_id"_2.fastq"  $sra_id"_R2".fastq
done

for file in *.fastq ; do
    echo "G-zipping $file" ;
    gzip $file ;
done