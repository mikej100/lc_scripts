#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=SRA-download-from-list
#SBATCH --ntasks=4
#SBATCH --mem=120G
#SBATCH --time=00-01:00:00

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
echo "starting SRA download script"

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

echo "SRA ids: "${sra_ids[@]}
for sra_id in ${sra_ids[@]}; do
    fastq-dump --split-files $sra_id
    mv $sra_id"_1.fast1  $model"_R1
    mv $sra_id"_2.fast1  $model"_R2
done

for file in *.fastq ; do
    echo "G-zipping $file" ;
    gzip $file ;
done


# for file in *.fastq.gz ; do echo "$file" ; done
