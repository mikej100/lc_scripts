#!/bin/bash
#SBATCH --partition=long
#SBATCH --ntasks=5
#SBATCH --mem=50G
#SBATCH --time=00-10:00:00
#SBATCH --output=slurm/%j_%x.out
#SBATCH --error=slurm/%j_%x.out
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
# Change the genome reference in the bedToBigBed to match your sample
#
### Run as the following replacing <Model> with your sample name (without .bam)

# sbatch Bigwigs_for_lanceotron.sh <Model> 

now() {
    date +"%Y-%m-%dT%T"
}
echo $(now) Starting $(basename "${BASH_SOURCE}")

# Show git info for scripts folder
${SCRIPTS}/scripts_info.sh || true
#Model name is the folder name
workingdir="$(pwd)"
model=$(echo $(pwd) | awk -F/ '{print $NF}') ;

## Generates a Bigwig of bin size 1 and RPKM normalised, the preferred input for Lanceotron
#bamCoverage --bam "$bam" -o "$model".bw --extendReads -bs 1 --normalizeUsing RPKM ;

bw_file="${model}.bw"
echo "$(now) Starting Lanceotron for ${bw_file}" ;

# Runs Lanceotron
lanceotron callPeaks $bw_file -f Lanceotron_output




L_tron_bed_fname="Lanceotron_output/${model}_L-tron.bed" 
L_tron_toppeaks_fname="Lanceotron_output/${model}_L-tron_toppeaks.bed"

awk 'BEGIN{FS=OFS="\t"} { if(($4>=0.8) && ($4<=1)) { print $1,$2,$3} }' \
 $L_tron_bed_fname > \
 $L_tron_toppeaks_fname

### Generating a bigbed for UCSC
L_tron_toppeaks_sorted_fname="Lanceotron_output/${model}_L-tron_toppeaks.sorted.bed"
module load ucsctools/385

echo "$(now) Starting sort for ${L_tron_toppeaks_fname}" ;
sort -k1,1 -k2,2n \
	$L_tron_toppeaks_fname > \
	$L_tron_toppeaks_sorted_fname
	
L_tron_toppeaks_bb_fname="Lanceotron_output/${model}_L-tron_toppeaks.bb"
echo "$(now) Starting bedToBigBed for ${L_tron_toppeaks_sorted_fname}"
bedToBigBed \
	$L_tron_toppeaks_sorted_fname \
 	$chr_sizes \
	$L_tron_toppeaks_bb_fname