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

#=== Boiler plate for job logging ==============================================
now() {
    date +"%Y-%m-%dT%T"
}
echo $(now) Starting $(basename "${BASH_SOURCE}")
# Show git info for scripts folder
${SCRIPTS}/scripts_info.sh || true
#===============================================================================

#endedness="paired"
endedness="single" 
genome="${1:-mm39}"
echo $genome



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
workingdir="$(pwd)"
model=$(echo $(pwd) | awk -F/ '{print $NF}') ;

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

# Read configuration data
mapfile -t sra_ids < sra_ids.cfg
mapfile -t repl < replicate_names.cfg
echo "SRA ids: "${sra_ids[@]}
echo "Replicate names: "${repl[@]}
echo "ref_genome: ${ref_genome}"

# Set extended read length for single-end reads
if [[ "${endedness}" == "single" ]]; then
    ext_readlength="200"
else    
    ext_readlength=""
fi

## Pipeline stesp ==============================================================
#        bowtie2 -k 2 -N 1 -p 4 -U  $repl \
#        bowtie2 -k 2 -N 1 -p 4  --sra-acc SRR5453542 \
#        bowtie2 -k 2 -N 1 -p 4   -1 ${repl}_R1.fastq.gz -2 ${repl}_R2.fastq.gz \

for r in ${repl[@]}; do
if false; then
    if [[ "${endedness}" == "single" ]]; then
        echo "$(now) Calling bowtie for ${repl} as single-ended data"
        echo "$(now) Calling bowtie for source file  ${sra_ids[0]}"
        bowtie2 -k 2 -N 1 -p 4 -U  ${sra_ids[0]} \
            --maxins 1000 \
            -x $ref_genome \
            -S ${r}_alignment.sam
    else
        echo "$(now) Calling bowtie for ${repl} as pair-ended data"
#        bowtie2 -k 2 -N 1 -p 4  --sra-acc ${repl}
        bowtie2 -k 2 -N 1 -p 4   -1 ${repl}_R1.fastq.gz -2 ${repl}_R2.fastq.gz \
            --maxins 1000 \
            -x $ref_genome \
            -S ${repl}_alignment.sam
    fi
        #Filter the mapped reads (filtering the mapped reads from unmapped) see https://www.htslib.org/doc/samtools-view.html for the explaination of flags.

    echo "$(now) Calling samtools::view for ${repl}"
        samtools view -@8 -bS -F4 -o ${repl}_mapped.bam ${repl}_alignment.sam ;

        #Sort the reads in mapped.bam by coordinate.
    echo "$(now) Calling samtools::sort for ${repl}"
        samtools sort ${repl}_mapped.bam -@8 -o ${repl}_mapped.srtC ;
        rm -rf ${repl}_mapped.bam

        #Remove PCR duplicates.
    echo "$(now) Calling samtools::rmdup for ${repl}"
        samtools rmdup ${repl}_mapped.srtC ${repl}_filtered.bam ; 
        rm -rf ${repl}_mapped.srtC *.txt *.histogram *.hist

        #Index the bam file
    echo "$(now) Calling samtools::index for ${repl}"
        samtools index ${repl}_filtered.bam  

    echo "$(now) renaming bam files for ${repl}"
        mv ${repl}_filtered.bam       ${repl}.bam
        mv ${repl}_filtered.bam.bai   ${repl}.bai

fi
    echo "$(now) Calling bamCoverage for ${repl}"
    bamCoverage --bam "${repl}.bam" \
            -o "${repl}.bw" \
            --extendReads ${ext_readlength} \
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
