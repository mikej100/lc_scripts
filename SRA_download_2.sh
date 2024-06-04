#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=SRA
#SBATCH --ntasks=4
#SBATCH --mem=120G
#SBATCH --time=00-23:59:59

module load sratoolkit
vdb-config --interactive

fastq-dump --split-files SRR15970891 ; #ter119_APH_spleen_rep1

mv SRR15970891_1.fastq SRR15970891_ter119_APH_spleen_rep1_R1.fastq ; SRR15970891_2.fastq SRR15970891_ter119_APH_spleen_rep1_R2.fastq ;
for file in *.fastq ; do
    echo "G-zipping $file" ;
    gzip $file ;
done


for file in *.fastq.gz ; do echo "$file" ; done
