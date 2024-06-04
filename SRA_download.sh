#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=SRA
#SBATCH --ntasks=4
#SBATCH --mem=120G
#SBATCH --time=06-23:59:59

module load sratoolkit

fastq-dump --split-files SRR12087137 ; #WT_rep1
fastq-dump --split-files SRR12087138 ; #WT_rep2
fastq-dump --split-files SRR12087139 ; #WT_rep3
fastq-dump --split-files SRR12087140 ; #D4448_rep1
fastq-dump --split-files SRR12087141 ; #D4448_rep2
fastq-dump --split-files SRR12087142 ; #D4448_rep3

mv SRR12087137_1.fastq SRR12087137_WT_rep1_R1.fastq ;
mv SRR12087137_2.fastq SRR12087137_WT_rep1_R2.fastq ;

mv SRR12087138_1.fastq SRR12087138_WT_rep2_R1.fastq ;
mv SRR12087138_2.fastq SRR12087138_WT_rep2_R2.fastq ;

mv SRR12087139_1.fastq SRR12087139_WT_rep3_R1.fastq ;
mv SRR12087139_2.fastq SRR12087139_WT_rep3_R2.fastq ;

mv SRR12087140_1.fastq SRR12087140_D4448_rep1_R1.fastq ;
mv SRR12087140_2.fastq SRR12087140_D4448_rep1_R2.fastq ;

mv SRR12087141_1.fastq SRR12087141_D4448_rep2_R1.fastq ;
mv SRR12087141_2.fastq SRR12087141_D4448_rep2_R2.fastq ;

mv SRR12087142_1.fastq SRR12087142_D4448_rep3_R1.fastq ;
mv SRR12087142_2.fastq SRR12087142_D4448_rep3_R2.fastq ;

for file in *.fastq ; do
    echo "G-zipping $file" ;
    gzip $file ;
done


for file in *.fastq.gz ; do echo "$file" ; done