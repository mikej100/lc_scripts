#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=SRA
#SBATCH --ntasks=4
#SBATCH --mem=120G
#SBATCH --time=00-00:01:00

sra_ids = '"SRR14573546"'
# sra_ids = '"SRR14573546"  "SRR14573547"  "SRR14573548"'

echo "Show git commit reference:"
git -C ./lc_scripts/ show -s --format=%h%x09%ci%x09%an%x09%s

#Model name is the folder name
model=$(echo $(pwd) | awk -F/ '{print $NF}') ;

#Create a date and timestamp for when analysis began
Start_time=`date`
Time=`date +"%T"`

module load sratoolkit
vdb-config --interactive

for sra_id in $sra_ids; do
    fastq-dump --split-files $sra_id
    mv $sra_id"_1.fast1  $model"_R1
    mv $sra_id"_2.fast1  $model"_R2
done

for file in *.fastq ; do
    echo "G-zipping $file" ;
    gzip $file ;
done


for file in *.fastq.gz ; do echo "$file" ; done
