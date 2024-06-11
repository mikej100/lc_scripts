#!/bin/bash
#SBATCH --partition=long
#SBATCH --ntasks=5

# This is a script to average 2 or 3 bigwigs e.g. two reps of the same ChIP or ATAC.
# Bigwigs must have already been generated
# Use in an environment with wiggletools and ucsc-wigtobigwig 
# mamba activate bed_sam_deep_tools

source /project/higgslab/lcornell/mamba_installation/conda/bin/activate bed_sam_deep_tools

## This is how to run: 
#sbatch bigwig_averaging_wiggle.sh <output file name _averaged> <Input 1.bw> <Input 2.bw> (optional) <Input 3.bw>
# Make sure the genome build is correct

# sbatch bigwig_averaging_wiggle.sh GSM923582_mm9_Tal1_ter119_averaged.bw GSM923582_mm9_wgEncodePsuTfbsErythroblTal1BE14halfCd1InputRepSignalRep1.bigWig GSM923582_mm9_wgEncodePsuTfbsErythroblTal1BE14halfCd1InputRepSignalRep2.bigWig GSM923582_mm9_wgEncodePsuTfbsErythroblTal1BE14halfCd1InputRepSignalRep3.bigWig
# sbatch bigwig_averaging_wiggle.sh APHSpleen_Ter119_ATAC_averaged.bs1.rpkm.bw APHSpleen_Ter119_ATAC_rep1.bs1.rpkm.bw APHSpleen_Ter119_ATAC_rep2.bs1.rpkm.bw

output_file="$1"
file_1="$2"
file_2="$3"
file_3="$4"

echo "Input bigwigs : "$file_1" "$file_2" "$file_3" " ;
echo 

echo "Starting Wiggletools mean"
echo

if [ -z $4 ]; then
	wiggletools mean "$file_1" "$file_2" > "$output_file".temp.wig
else
	wiggletools mean "$file_1" "$file_2" "$file_3" > "$output_file".temp.wig
fi ;

echo "Starting wigToBigWig"
echo 

wigToBigWig "$output_file".temp.wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt "$output_file"

rm "$output_file".temp.wig

echo "Finished"
echo 

exit ;
 
