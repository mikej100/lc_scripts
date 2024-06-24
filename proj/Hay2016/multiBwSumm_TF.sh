#!/bin/bash
#SBATCH --partition=short
#SBATCH --output=slurm/%j_%x.out 
#SBATCH --error=slurm/%j_%x.out

###############################################################################
###############################################################################
descr=\
"MultiBigWig summary for tranlsation factors: 500bp slop added to ATAC bed"
#   Run this from an analysis directory at same level as assay directories.
###############################################################################
#                 In-code configuration
slop=500
bed_filename="../ATAC/Lanceotron_output/ATAC_L-tron_toppeaks_slop${slop}.bed"
out_file="mbws_cf"
###############################################################################
#                         Boiler plate for job logging
now() {
    date +"%Y-%m-%dT%T"
}
echo $(now) Starting $(basename "${BASH_SOURCE}")
# Show git info for scripts folder
${SCRIPTS}/scripts_info.sh || true
###############################################################################
###############################################################################
#                         Set up configuration
#Model name is the folder name
echo ${descr}
echo "Bigwig filenames: ${bw_filenames[@]}"

echo "bed filename ${bed_filename}"

echo "Assay names: ${assay[@]}"
echo "labels: ${label[@]}"

module purge
module load python-cbrg/current 
#===============================================================================
#                  Run process steps
#
echo "$(now) Starting multiBigwigSummary"
multiBigwigSummary BED-file \
     -b \
        "../Med1/Med1.bw" \
     -o ${out_file}.npz \
     --outRawCounts ${out_file}.tab \
     --labels Med_1 \
     --BED $bed_filename
# 