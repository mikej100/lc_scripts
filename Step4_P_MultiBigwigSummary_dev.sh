#!/bin/bash
#SBATCH --partition=short
#SBATCH --output=slurm/%j_%x.out 
#SBATCH --error=slurm/%j_%x.out

# Instructions to run
# Create directory for he analysis at the same level as the models
# cd to the analysis directory and run this script
#=== Boiler plate for job logging ==============================================
now() {
    date +"%Y-%m-%dT%T"
}
echo $(now) Starting $(basename "${BASH_SOURCE}")
# Show git info for scripts folder
${SCRIPTS}/scripts_info.sh || true
#===============================================================================

#==================== Configuration section ====================================
# Get the list of models

# Get the models to include in the analysis
# TODO if file not found then use the models.cfg in the project directory
mapfile -t model < models.cfg
models=${model[@]}
echo "Models: ${model[@]}"

bw_filenames=()
for i in ${!model[@]}; do
     bw_filenames+="../${model[i]}/${model[i]}.bw "
done
echo "Bigwig filenames: ${bw_filenames}"

bed_filename="../${model[0]}/Lanceotron_output/${model[0]}_L-tron_toppeaks.bed"
echo "bed filename ${bed_filename}"
 
#===============================================================================
echo "$(now) Starting multiBigwigSummary"
multiBigwigSummary BED-file \
     -b $bw_filenames \
     -o mbw_results.npz \
     --outRawCounts mbw_results.tab \
     -l $models \
     --BED $bed_filename
# 
# 
# 
# 
# #  multiBigwigSummary BED-file \
# #  -b /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/H3K4me1.bw \
# #  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/H3K4me3.bw \
# #  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Gata1.bw \
# #  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Klf1.bw \
# #  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Nfe2.bw \
# #  /project/higgslab/lcornell/Super_enhancer_analysis/Bigwigs/Med1.bw \
# #  -o enhancer_mcs_mm9.npz --outRawCounts enhancer_mcs_mm9.tab \
# #  --labels H3K4me1 H3K4me3 Gata1 Klf1 Nfe2 Med1 \
# #  --BED /project/higgslab/lcornell/BEDs/SuperEnhancers/enhancer_mcs_mm9.bed
# # 