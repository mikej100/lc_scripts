#!/bin/bash
#SBATCH --partition=long
###SBATCH --job-name=SRA
#SBATCH --ntasks=5
#SBATCH --mem=10G
#SBATCH --time=00-30:00:00
#SBATCH --output=slurm/%j_%x.out 
#SBATCH --error=slurm/%j_%x.out
###SBATCH --mail-user=
##SBATCH --mail-type=fail
# New line by Mike 20240619T1705
# Based on ChiiP_ATAC_pipeline_Bowtie2_fast.sh script from Lucy Cornell, modified by Mike Jennings for use with
# lanceotron.
# Changes are:
#       not having 8x lane fastq from same samples. 
#          I need to understand this requirement
#       skipped the bamCoverage, as this step is in next script.
####### If you want to mask blacklist regions, see line 190 

####### To RUN: sbatch ChiP_ATAC_pipeline_Bowtie2_fast.sh <genome>

###############################################################################
#                         Get command line options
# Function to display help text
usage() {
    echo "Usage: $(basename $0) [-r N,...] [-U] "
    echo "Options:"
    echo "  -g          Genome alignment. Default is mm39"
    echo "  -h          Display this help message"
    echo "  -r          select replicates to process, by postion. "
    echo "              comma delimited no spaces. E.g. -r 2,3 "
    echo "              Default is all items in config file."
    echo "  -U          Unpaired read (single-ended). Default is paired."
    echo "  -V          very verbose: print every command line"
}
# Defaults
endedness="paired"
genome="mm39"

# Parse options using getopts
while getopts "hr:Uv" option; do
    case "${option}" in
        g)  genome=($OPTARG)
            ;;
        r)  IFS=,
            repl_indices=($OPTARG)
            ;;
        U)  endedness="single"
            ;;
        h)  # Help option
            usage
            exit 0
            ;;
        V)  set -x
            ;;
        \?) # Invalid option
            echo "Invalid option: -$OPTARG"
            usage
            exit 1
            ;;
    esac
done
shift $(($OPTIND - 1))

###############################################################################
echo "Options settings:"
echo "  genome: ${genome}"
echo "  endedness: ${endedness}"
echo "  replicate indices: ${repl_indices[@]}"
echo "  Remaining arguments read: $*"
echo 
###############################################################################
###############################################################################
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
model="$(basename $PWD)"

# Read configuration data
mapfile -t repl < replicate_names.cfg
echo "Replicate names: "${repl[@]}

#Generally, genome will be genome will be mm10 or hg19. Can use if statement to set correct species based on genome build.
if [ "$genome" == "mm9" ] || [ "$genome" == "mm10" ] || [ "$genome" == "mm39" ]; then 
        species="Mus_musculus"
elif [ "$genome" == "hg19" ] || [ "$genome" == "hg38" ] ; then 
        species="Homo_sapiens"
else
	echo "Please check your genome build and/or species."  
	exit
fi
ref_genome="/databank/igenomes/"$species"/UCSC/"$genome"/Sequence/Bowtie2Index/genome"
echo "ref_genome: ${ref_genome}"

module purge #Unloads all modules from the users environment
#Initialise the required modules/packages for analysis.
module load samtools/1.17 bowtie2/2.4.2 bedtools/2.25.0 ucsctools/385 trim_galore/0.6.10 flash/1.2.11;
module list #Prints list of modules and versions


################################################################################
#                       Run process steps
#                                                                              #
#Use bowtie2 to align the fastq files to the designated genome (set as $genome).
#Currently, allow reads which multi-map to max 2 regions to continue into downstream analysis. Change -k to alter this.
#Currently, allow maximum of one mismatch in seed. Change -N to alter this. Allowing fewer mismatches will increase accuracy; allowing more
#will increase sensitivity.
#Currently, maximum fragment length allowed is 1000bp. Change --maxins if you want to alter this e.g. perhaps change to 150bp if you want
#only mononucleosomes.

# Use extendRead option only for paired ended data.
if [[ "${endedness}" == "single" ]]; then
    extendRead_option=""
else    
    extendRread_option="--extendReads "
fi

# if [ -z $repl_indices]; do
#     for i in ${!repl_indices[@]}: do
#         repl_todo+=$repl[$i]
#     done;
# else
#     repl_todo=repl
# fi
# 
# echo "${repl_todo[']}"
# exit
## Pipeline stesp ==============================================================
#        bowtie2 -k 2 -N 1 -p 4 -U  $repl \
#        bowtie2 -k 2 -N 1 -p 4  --sra-acc SRR5453542 \
#        bowtie2 -k 2 -N 1 -p 4   -1 ${_repl}_R1.fastq.gz -2 ${_repl}_R2.fastq.gz \

for _repl in ${repl[@]}; do
    if [[ "${endedness}" == "single" ]]; then
        echo "$(now) Calling bowtie for ${_repl} as single-ended data"
        bowtie2 -k 2 -N 1 -p 4 -U  ${_repl} \
            --maxins 1000 \
            -x $ref_genome \
            -S ${r}_alignment.sam
    else
        echo "$(now) Calling bowtie for ${_repl} as pair-ended data"
#        bowtie2 -k 2 -N 1 -p 4  --sra-acc ${_repl}
        bowtie2 -k 2 -N 1 -p 4   -1 ${_repl}_R1.fastq.gz -2 ${_repl}_R2.fastq.gz \
            --maxins 1000 \
            -x $ref_genome \
            -S ${_repl}_alignment.sam
    fi
        #Filter the mapped reads (filtering the mapped reads from unmapped) see https://www.htslib.org/doc/samtools-view.html for the explaination of flags.

    echo "$(now) Calling samtools::view for ${_repl}"
        samtools view -@8 -bS -F4 -o ${_repl}_mapped.bam ${_repl}_alignment.sam ;

        #Sort the reads in mapped.bam by coordinate.
    echo "$(now) Calling samtools::sort for ${_repl}"
        samtools sort ${_repl}_mapped.bam -@8 -o ${_repl}_mapped.srtC ;
        rm -rf ${_repl}_mapped.bam

        #Remove PCR duplicates.
    echo "$(now) Calling samtools::rmdup for ${_repl}"
        samtools rmdup ${_repl}_mapped.srtC ${_repl}_filtered.bam ; 
        rm -rf ${_repl}_mapped.srtC *.txt *.histogram *.hist

        #Index the bam file
    echo "$(now) Calling samtools::index for ${_repl}"
        samtools index ${_repl}_filtered.bam  

    echo "$(now) renaming bam files for ${_repl}"
        mv ${_repl}_filtered.bam       ${_repl}.bam
        mv ${_repl}_filtered.bam.bai   ${_repl}.bai

    echo "$(now) Calling bamCoverage for ${_repl}"
    bamCoverage --bam "${_repl}.bam" \
            -o "${_repl}.bw" \
            ${extend_read_option} \
            -bs 1 \
            --normalizeUsing RPKM
done

module unload python-base
module load python-cbrg

#Create the bigwig (CPM (can swap CPM for RPKM here if you want. Other normalisation options available too) and smoothed over 100bp in
#50bp bins. 50bp bins is default but can change by altering -bs flag. Can change smoothing window by changing --smoothLength. Can also change
#normalisation method by altering --normalizeUsing. Increase -bs and --smoothLength metrics for smoother data (less noisy),
#Decrease both to improve resolution...better if sequencing is deeper

blacklist=$(echo $(ls /project/higgslab/shared/00_analysis_file_bank/Blacklists/$genome*.bed))

# If you would like blacklist regions removed, add this to the bamcoverage command:
# --blackListFileName "$blacklist"
# # Change names
#Get rid of intermediate files to save space. But KEEP filtered.bam for Macs2 analysis (need bam for peak calling and downstream analysis).
rm -rf  *_filtered.bg \
        *_filtered.srt.bg \
        *_mapped.srtC.bam \
        *_alignment.sam ;
rm *fq*

echo -e "$(now) Intermediate files removed at `date +"%T"`\n" ;
echo
echo -e "$(now) ChIP pipeline finished `date +"%T"`"
echo
exit ;
