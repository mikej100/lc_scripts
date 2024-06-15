#!/bin/bash
#SBATCH --partition=long
###SBATCH --job-name=SRA
#SBATCH --ntasks=5
#SBATCH --mem=10G
#SBATCH --time=00-10:00:00
#SBATCH --output=slurm/%j_%x.out 
#SBATCH --error=slurm/%j_%x.out
###SBATCH --mail-user=
##SBATCH --mail-type=fail
# Based on ChiiP_ATAC_pipeline_Bowtie2_fast.sh script from Lucy Cornell, modified by Mike Jennings for use with
# lanceotron.
# Changes are:
#       not having 8x lane fastq from same samples. 
#          I need to understand this requirement
#       skipped the bamCoverage, as this step is in next script.
####### If you want to mask blacklist regions, see line 190 

####### To RUN: sbatch ChiP_ATAC_pipeline_Bowtie2_fast.sh <genome>
###########################################################################################################################################################
#move to working directory

####### Set Opions 
###############################################################################
#endedness="paired"
endedness="single" 

genome="${1:-mm39}"
echo $genome


#Create a date and timestamp for when analysis began
Start_time=`date`
Time=`date +"%T"`
workingdir="$(pwd)"


now() {
    date +"%Y-%m-%dT%T"
}


echo $(now) Starting $(basename "${BASH_SOURCE}")
# Show git info for scripts folder
${SCRIPTS}/scripts_info.sh || true
# echo
# echo "Script file information from git"
# echo "================================"
# 
# echo "Here are repository; branch; most recent commit reference, date, author and comment"
# 
# SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# git -C "$SCRIPT_DIR" remote -v | grep fetch
# # git -C $(dirname ${BASH_SOURCE}) rev-parse --abbrev-ref HEAD
# # git -C $(dirname ${BASH_SOURCE}) show -s --format=%h%x09%ci%x09%an%x09%sk


#git -C "${SCRIPT_DIR}ASH_SOURCE})" rev-parse --abbrev-ref HEAD  
#Genome path, default to mm39

#Generally, genome will be genome will be mm10 or hg19. Can use if statement to set correct species based on genome build.
if [ "$genome" == "mm9" ] || [ "$genome" == "mm10" ] || [ "$genome" == "mm39" ]; then 
        species="Mus_musculus"
elif [ "$genome" == "hg19" ] || [ "$genome" == "hg38" ] ; then 
        species="Homo_sapiens"
else
	echo "Please check your genome build and/or species."  
	exit
fi ;
ref_genome="/databank/igenomes/"$species"/UCSC/"$genome"/Sequence/Bowtie2Index/genome"

#Model name is the folder name
model=$(echo $(pwd) | awk -F/ '{print $NF}') ;

#move to directory

echo "----------------- Bowtie2 Alignment for ChIP ----------------- " ;
echo
echo "Analysis started: $Start_time" ;
echo
echo "Current working directory is: "$workingdir" "
echo
echo "Genome Bowtie2 index: "$species" "$genome" "
echo
echo "Sample: "${model}" "
echo

module purge #Unloads all modules from the users environment
#Initialise the required modules/packages for analysis.
module load samtools/1.17 bowtie2/2.4.2 bedtools/2.25.0 ucsctools/385 trim_galore/0.6.10 flash/1.2.11;
module list #Prints list of modules and versions



#Use bowtie2 to align the fastq files to the designated genome (set as $genome).
#Currently, allow reads which multi-map to max 2 regions to continue into downstream analysis. Change -k to alter this.
#Currently, allow maximum of one mismatch in seed. Change -N to alter this. Allowing fewer mismatches will increase accuracy; allowing more
#will increase sensitivity.
#Currently, maximum fragment length allowed is 1000bp. Change --maxins if you want to alter this e.g. perhaps change to 150bp if you want
#only mononucleosomes.

# Read SRA_ids from a a config file. One per line.
mapfile -t sra_ids < sra_ids.cfg
echo "SRA ids: "${sra_ids[@]}

echo "ref_genome: ${ref_genome}"
for sra_id in ${sra_ids[@]}; do
    if [[ "${endedness}" == "single" ]]; then
        echo "$(now) Calling bowtie for ${sra_id} as single-ended data"
        bowtie2 -k 2 -N 1 -p 4 -U  --sra-acc $sra_id \
#        bowtie2 -k 2 -N 1 -p 4 -U  $sra_id \
            --maxins 1000 \
            -x $ref_genome \
            -S ${outfile}_alignment.sam
    else
        echo "$(now) Calling bowtie for ${sra_id} as pair-ended data"
#        bowtie2 -k 2 -N 1 -p 4  --sra-acc ${sra_id}
        bowtie2 -k 2 -N 1 -p 4  c ${sra_id}
            --maxins 1000 \
            -x $ref_genome \
            -S ${sra_id}_alignment.sam
    fi
        #Filter the mapped reads (filtering the mapped reads from unmapped) see https://www.htslib.org/doc/samtools-view.html for the explaination of flags.

    echo "$(now) Calling samtools::view for ${sra_id}"
        samtools view -@8 -bS -F4 -o ${sra_id}_mapped.bam ${sra_id}_alignment.sam ;

        #Sort the reads in mapped.bam by coordinate.
    echo "$(now) Calling samtools::sort for ${sra_id}"
        samtools sort ${sra_id}_mapped.bam -@8 -o ${sra_id}_mapped.srtC ;
        rm -rf ${sra_id}_mapped.bam

        #Remove PCR duplicates.
    echo "$(now) Calling samtools::rmdup for ${sra_id}"
        samtools rmdup ${sra_id}_mapped.srtC ${sra_id}_filtered.bam ; 
        rm -rf ${sra_id}_mapped.srtC *.txt *.histogram *.hist

        #Index the bam file
    echo "$(now) Calling samtools::index for ${sra_id}"
        samtools index ${sra_id}_filtered.bam  

    echo "$(now) renaming bam files for ${sra_id}"
        mv ${sra_id}_filtered.bam       ${sra_id}.bam
        mv ${sra_id}_filtered.bam.bai   ${sra_id}.bai
        mv ${sra_id}_filtered.rpkm.bw   ${sra_id}.rpkm.bw

    echo "$(now) Calling bamCoverage for ${sra_id}"
    
    bamCoverage --bam "${sra_id}.bam" \
            -o "${sra_id}.bw" \
            --extendReads \
            -bs 1 \
            --normalizeUsing RPKM
done

module unload python-base

echo
echo "----------------- Bamcoverage for Bigwigs -----------------" ;
echo

module load python-cbrg

#Create the bigwig (CPM (can swap CPM for RPKM here if you want. Other normalisation options available too) and smoothed over 100bp in
#50bp bins. 50bp bins is default but can change by altering -bs flag. Can change smoothing window by changing --smoothLength. Can also change
#normalisation method by altering --normalizeUsing. Increase -bs and --smoothLength metrics for smoother data (less noisy),
#Decrease both to improve resolution...better if sequencing is deeper

blacklist=$(echo $(ls /project/higgslab/shared/00_analysis_file_bank/Blacklists/$genome*.bed))

# If you would like blacklist regions removed, add this to the bamcoverage command:
# --blackListFileName "$blacklist"
# Mike Jennings comment out this block becuase the .bw files are generated in "Step1" script"
# bamCoverage -b filtered.bam --normalizeUsing RPKM --smoothLength 100 -p 4 -bs 50 -o filtered.rpkm.bw ;
# 
# if [ -s filtered.rpkm.bw ]
#         then
#         echo "RPKM normalized bigwig successfully generated at `date +"%T"`." 
#         echo "See results in filtered.rpkm.bw."
#         echo
#         else
#         echo "bigwig generation unsuccessful."  
#         exit
# fi
# 
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
