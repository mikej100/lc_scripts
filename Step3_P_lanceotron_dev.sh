#!/bin/bash
#SBATCH --partition=long
#SBATCH --ntasks=5
#SBATCH --mem=50G
#SBATCH --time=00-10:00:00
#SBATCH --output=slurm/%j_%x.out
#SBATCH --error=slurm/j_%x.out
###SBATCH --mail-type=end,fail

chr_sizes_mm9="/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt"
chr_sizes_mm39="/databank/igenomes/Mus_musculus/UCSC/mm39/Sequence/WholeGenomeFasta/chr_sizes.txt"
chr_sizes=$chr_sizes_mm39

module purge
module load lanceotron
# Generate 
## This is an example of how to generate bigwigs and peak call fromt hese peaks using lanceotron
# This version by Mike Jennings based on script from Lucy Cornell
# Changed to run from within the model folder
# 
## This is also written for samples which are held in individual folders, that have files of the same name
## ie 
# ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2.bam
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2.bam.bai
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2_R1.fastq.gz
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2_R2.fastq.gz
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2.smooth.rpkm.bw
#
## This should also be run from within the sample directory ie 
# ├ Experiment folder					<------ From here
# ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2.bam
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2.bam.bai
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2_R1.fastq.gz
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2_R2.fastq.gz
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2.smooth.rpkm.bw
#
# Change the genome reference in the bedToBigBed to match your sample
#
### Run as the following replacing <Model> with your sample name (without .bam)

# sbatch Bigwigs_for_lanceotron.sh <Model> 
#Create a date and timestamp for when analysis began

now() {
    date +"%Y-%m-%dT%T"
}
NOW=`date +"%Y-%m-%dT%T"`

echo $(now) Starting $(basename "${BASH_SOURCE}")

workingdir="$(pwd)"

# Show git info for scripts folder
${SCRIPTS}/scripts_info.sh || true
#Model name is the folder name
model=$(echo $(pwd) | awk -F/ '{print $NF}') ;

# Find the bam with the corresponding name within its folder
# 
#echo "Starting Bigwig for $bam" ;
#
## Generates a Bigwig of bin size 1 and RPKM normalised, the preferred input for Lanceotron
#bamCoverage --bam "$bam" -o "$model".bw --extendReads -bs 1 --normalizeUsing RPKM ;

input_file="${model}.bw"
echo "$(now) Starting Lanceotron for ${input_file}" ;

# Runs Lanceotron
lanceotron callPeaks $input_file -f Lanceotron_output

# This filters peak calls to only those of a score greater than 0.8
echo "Taking top scoring peaks" ;

awk 'BEGIN{FS=OFS="\t"} { if(($4>=0.8) && ($4<=1)) { print $1,$2,$3,$4,$5,$6,$7,$8 } }' Lanceotron_output/"$model"_L-tron.bed > Lanceotron_output/"$model"_L-tron_toppeaks.bed ;

### Generating a bigbed for UCSC
module load ucsctools/385

bedToBigBed Lanceotron_output/"$model"_L-tron_toppeaks.bed \
 	$chr_sizes \
	Lanceotron_output/"$model"_L-tron_toppeaks_bigGenePred.bb ;