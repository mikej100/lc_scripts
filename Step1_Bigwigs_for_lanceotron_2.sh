#!/bin/bash
#SBATCH --partition=long
#SBATCH --ntasks=5
#SBATCH --mem=50G
#SBATCH --time=00-10:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err
#SBATCH --mail-type=end,fail

module load python-cbrg/202401
module load lanceotron/20230726

## This is an example of how to generate bigwigs and peak call fromt hese peaks using lanceotron
## This is also written for samples which are held in individual folders, that have files of the same name
## ie 
# ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2.bam
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2.bam.bai
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2_R1.fastq.gz
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2_R2.fastq.gz
# │   ├── CD71_EB_Rad21_Mira_INVEN_E14_Rep2.smooth.rpkm.bw
#
## This should also be run from outside the sample directory ie 
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
Start_time=`date`
Time=`date +"%T"`
workingdir="$(pwd)"
echo "Show git commit reference:"
git -C ./lc_scripts/ show -s --format=%h%x09%ci%x09%an%x09%s

# Model name from the first input into the script
if ! model="$1"; then
	echo "Add model name"
	exit 1
fi

# Find the bam with the corresponding name within its folder
bam="$model"/"$model".bam
# 
echo "Starting Bigwig for $bam" ;

# Generates a Bigwig of bin size 1 and RPKM normalised, the preferred input for Lanceotron
bamCoverage --bam "$bam" -o Peak_analysis/"$model".bw --extendReads -bs 1 --normalizeUsing RPKM ;

echo "Starting Lanceotron for $bam" ;

# Runs Lanceotron
lanceotron callPeaks Peak_analysis/"$model".bw -f Peak_analysis/Lanceotron_outputs ;

# This filters peak calls to only those of a score greater than 0.8
echo "Taking top scoring peaks" ;

awk 'BEGIN{FS=OFS="\t"} { if(($4>=0.8) && ($4<=1)) { print $1,$2,$3,$4,$5,$6,$7,$8 } }' Peak_analysis/Lanceotron_outputs/"$model"_L-tron.bed > Peak_analysis/Lanceotron_outputs/"$model"_L-tron_toppeaks.bed ;

### Generating a bigbed for UCSC
module load ucsctools/385

bedToBigBed Peak_analysis/Lanceotron_outputs/"$model"_L-tron_toppeaks.bed /databank/igenomes/Mus_musculus/UCSC/mm39/Sequence/WholeGenomeFasta/chr_sizes.txt Peak_analysis/Lanceotron_outputs/"$model"_L-tron_toppeaks_bigGenePred.bb ;

