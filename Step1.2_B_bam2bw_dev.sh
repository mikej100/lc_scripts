#!/bin/bash
#SBATCH --partition=long
###SBATCH --job-name=SRA
#SBATCH --ntasks=5
#SBATCH --mem=10G
#SBATCH --time=00-10:00:00
#SBATCH --output=slurm/%j_%x.out ### If does not exist job will fail with exit code 53
### Send err to same file for debugging use
#SBATCH --error=slurm/%j_%x.out
###SBATCH --mail-user=
#SBATCH --mail-type=fail

# subset of step 1, to process mini-bam file so bw for the mini project
source /project/higgslab/lcornell/mamba_installation/conda/bin/activate bed_sam_deep_tools

NOW=`date +"%Y-%m-%dT%T"`

echo ${NOW} Starting $(basename "${BASH_SOURCE}")
echo  " "
${SCRIPTS}/scripts_info.sh || true

# Read SRA_ids from a a config file. One per line.
mapfile -t sra_ids < sra_ids.cfg
echo "SRA ids: "${sra_ids[@]}

for sra_id in ${sra_ids[@]}; do
        #Index the bam file
    echo "${NOW} Calling samtools::index for ${sra_id}"
        samtools index ${sra_id}.bam  

   #  echo "${NOW} Calling bamCoverage for ${sra_id}"
   #  
   #  bamCoverage --bam "${sra_id}.bam" \
   #          -o "${sra_id}.bw" \
   #          --extendReads \
   #          -bs 1 \
   #          --normalizeUsing RPKM
done