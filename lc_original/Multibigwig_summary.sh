#!/bin/bash
#SBATCH --partition=short

module load python-cbrg/current

multiBigwigSummary BED-file \
 -b mm10.60way.phyloP60way.bw \
 -o enhancer_mcs_mm10_conservation.npz --outRawCounts enhancer_mcs_mm10_conservation.tab \
 --labels phyloP --BED /project/higgslab/lcornell/BEDs/SuperEnhancers/enhancer_mcs_mm10.bed

 multiBigwigSummary BED-file \
 -b /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/H3K4me1.bw \
 /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/H3K4me3.bw \
 /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Gata1.bw \
 /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Klf1.bw \
 /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Nfe2.bw \
 /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Med1.bw \
 -o enhancer_mcs_mm9.npz --outRawCounts enhancer_mcs_mm9.tab \
 --labels H3K4me1 H3K4me3 Gata1 Klf1 Nfe2 Med1 \
 --BED /project/higgslab/lcornell/BEDs/SuperEnhancers/enhancer_mcs_mm9.bed