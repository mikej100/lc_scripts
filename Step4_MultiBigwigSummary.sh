#!/bin/bash
#SBATCH --partition=short

echo "Starting multiBigwigSummary"

module load python-cbrg/current
source /project/higgslab/lcornell/mamba_installation/conda/bin/activate bed_sam_deep_tools

model=("ATAC" "H3K4me1" "H3K4me3" "Med1")

echo $model

for i in "${!model[@]}"; do
    echo ${model[$i]}
    bw_f[$i]="../"${model[$i]}"/"${model[$i]}".bw"
done

for f in "${bw_f[@]}"; do
    echo $f
done



#  multiBigwigSummary BED-file \
#  -b /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/H3K4me1.bw \
#  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/H3K4me3.bw \
#  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Gata1.bw \
#  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Klf1.bw \
#  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Nfe2.bw \
#  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Med1.bw \
#  -o enhancer_mcs_mm9.npz --outRawCounts enhancer_mcs_mm9.tab \
#  --labels H3K4me1 H3K4me3 Gata1 Klf1 Nfe2 Med1 \
#  --BED /project/higgslab/lcornell/BEDs/SuperEnhancers/enhancer_mcs_mm9.bed
# 