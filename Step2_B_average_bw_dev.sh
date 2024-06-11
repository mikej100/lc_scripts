#!/bin/bash
#SBATCH --partition=long
#SBATCH --ntasks=5
#SBATCH --output=slurm/%j_%x.out 
#SBATCH --error=slurm/%j_%x.out

# This is a script to average 2 or 3 bigwigs e.g. two reps of the same ChIP or ATAC.
# Bigwigs must have already been generated
# Use in an environment with wiggletools and ucsc-wigtobigwig 
# mamba activate bed_sam_deep_tools


chr_sizes_mm9="/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt"
chr_sizes_mm39="/databank/igenomes/Mus_musculus/UCSC/mm39/Sequence/WholeGenomeFasta/chr_sizes.txt"
chr_sizes=$chr_sizes_mm39

source /project/higgslab/lcornell/mamba_installation/conda/bin/activate bed_sam_deep_tools

now() {
    date +"%Y-%m-%dT%T"
}

echo ${now} Starting $(basename "${BASH_SOURCE}")
# Show git info for scripts folder
#${SCRIPTS}/scripts_info.sh || true

echo "${now} DEBUG returned to main script"


#Model name is the folder name
model=$(echo $(pwd) | awk -F/ '{print $NF}') ;

## This is how to run: 
#sbatch bigwig_averaging_wiggle.sh <output file name _averaged> <Input 1.bw> <Input 2.bw> (optional) <Input 3.bw>
# Make sure the genome build is correct

# sbatch bigwig_averaging_wiggle.sh GSM923582_mm9_Tal1_ter119_averaged.bw GSM923582_mm9_wgEncodePsuTfbsErythroblTal1BE14halfCd1InputRepSignalRep1.bigWig GSM923582_mm9_wgEncodePsuTfbsErythroblTal1BE14halfCd1InputRepSignalRep2.bigWig GSM923582_mm9_wgEncodePsuTfbsErythroblTal1BE14halfCd1InputRepSignalRep3.bigWig
# sbatch bigwig_averaging_wiggle.sh APHSpleen_Ter119_ATAC_averaged.bs1.rpkm.bw APHSpleen_Ter119_ATAC_rep1.bs1.rpkm.bw APHSpleen_Ter119_ATAC_rep2.bs1.rpkm.bw


# Read SRA_ids from a a config file. One per line.
mapfile -t sra_ids < sra_ids.cfg
echo "SRA ids: "${sra_ids[@]}
output_file="$model"

fnames=$(printf %s.bw" " "${sra_ids[@]}")

# 
# file_1="$2"
# file_2="$3"
# file_3="$4"
# 
echo "Input bigwigs : $fnames"
# 
echo "${now} Starting Wiggletools mean"
wiggletools mean ${fnames} > "$output_file".temp.wig 
# 	wiggletools mean "$file_1" "$file_2" > "$output_file".temp.wig
# else
# 	wiggletools mean "$file_1" "$file_2" "$file_3" > "$output_file".temp.wig
# fi ;
# 
echo "${now} Starting wigToBigWig"
wigToBigWig "$output_file".temp.wig $chr_sizes "$output_file" ]
# 
# rm "$output_file".temp.wig
# 
 echo "${now}Finished"
 echo 