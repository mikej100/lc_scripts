#!/bin/bash
#SBATCH --partition=long
#SBATCH --ntasks=5
#SBATCH --mem=50G
#SBATCH --time=00-10:00:00
###SBATCH --output=%j_%x.out
###SBATCH --error=%j_%x.err
###SBATCH --mail-user=
#SBATCH --mail-type=end,fail

# Based on script from Lucy Cornell, modified by Mike Jennings for use with
# lanceotron.
# Changes are:
#       not having 8x lane fastq from same samples. 
#          I need to understand this requirement
#       skipped the bamCoverage, as this step is in next script.

## NEEDS TO HAVE 8x LANE fastq.gz from the same sample from FILES IN SAME FOLDER TO RUN ###
## 4 sequencing lanes for paired end reads

### Things to check before run:
####### Change working directory to where the fastqs are, please note the name of the directory is used to identify output files so make this name informative
####### Provide the genome build the analysis is running on
####### Check the module versions, early versions of modules may not work with this script
####### If you want to mask blacklist regions, see line 190 

####### To RUN: sbatch ChiP_ATAC_pipeline_Bowtie2_fast.sh <genome>
###########################################################################################################################################################
#move to working directory

#Create a date and timestamp for when analysis began
Start_time=`date`
Time=`date +"%T"`
workingdir="$(pwd)"

#Genome path, default to mm39
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

## Print inputs
echo -e "Input fastqs:\n$(ls *.fastq.gz | xargs -n1)" ;
echo

#Concatenate _R1 and _R2 fastq files into a single R1 and a single R2 file.
nohup cat *_R1* > R1.fastq.gz
nohup cat *_R2* > R2.fastq.gz

if [ -s R1.fastq.gz ]
	then
	echo "Read 1 reads successfully concatenated at `date +"%T"`."
        echo
	else
	echo "R1 concatenation unsuccessfull."  
	exit
fi

if [ -s R2.fastq.gz ]
        then
        echo "Read 2 reads successfully concatenated at `date +"%T"`."
        echo
        else
        echo "R2 concatenation unsuccessfull."  
        exit
fi

#Use bowtie2 to align the fastq files to the designated genome (set as $genome).
#Currently, allow reads which multi-map to max 2 regions to continue into downstream analysis. Change -k to alter this.
#Currently, allow maximum of one mismatch in seed. Change -N to alter this. Allowing fewer mismatches will increase accuracy; allowing more
#will increase sensitivity.
#Currently, maximum fragment length allowed is 1000bp. Change --maxins if you want to alter this e.g. perhaps change to 150bp if you want
#only mononucleosomes.

bowtie2 -k 2 -N 1 -p 4 \ -1 R1.fastq.gz \ -2 R2.fastq.gz \
--maxins 1000 \
-x /databank/igenomes/"$species"/UCSC/"$genome"/Sequence/Bowtie2Index/genome \
-S alignment.sam

if [ -s alignment.sam ]
        then
        echo "Alignment successful at `date +"%T"`." 
        echo "See results in alignment.sam."
        echo
        else
        echo "Alignment unsuccessful."  
        exit
fi

rm -rf *.fastq

#Filter the mapped reads (filtering the mapped reads from unmapped) see https://www.htslib.org/doc/samtools-view.html for the explaination of flags.

samtools view -@8 -bS -F4 -o mapped.bam alignment.sam ;

if [ -s mapped.bam ]
        then
        echo "Read filtering successful at `date +"%T"`.\n"
        echo "See results in mapped.bam."
        echo
        else
        echo "Read filtering unsuccessful."   
        exit
fi

#Sort the reads in mapped.bam by coordinate.

samtools sort mapped.bam -@8 -o mapped.srtC ;

if [ -s mapped.srtC ]
        then
        echo "Read sorting successful at `date +"%T"`. "
        echo "See results in mapped.srtC."
        echo
        else
        echo "Read sorting unsuccessful."  
        exit
fi

rm -rf mapped.bam

#Remove PCR duplicates.

samtools rmdup mapped.srtC filtered.bam ; 

if [ -s filtered.bam ]
        then
        echo "PCR duplicates removal successful at `date +"%T"`. "
        echo "See results in filtered.bam."
        echo
        else
        echo "PCR duplicates removal unsuccessful."   
        exit
fi

rm -rf mapped.srtC *.txt *.histogram *.hist

#Index the bam file
samtools index filtered.bam ; 

if [ -s filtered.bam.bai ]
        then
        echo "filtered.bam successfully indexed at `date +"%T"`. "
        echo "See results in filtered.bam.bai."
        echo
        else
        echo "bam indexing unsuccessful."   
        exit
fi

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
mv filtered.bam ${model}.bam
mv filtered.bam.bai ${model}.bam.bai
mv filtered.rpkm.bw ${model}.rpkm.bw

echo "Files renamed"
echo

#Get rid of intermediate files to save space. But KEEP filtered.bam for Macs2 analysis (need bam for peak calling and downstream analysis).
rm -rf  filtered.bg filtered.srt.bg mapped.srtC.bam alignment.sam ;
rm *fq*

echo -e "Intermediate files removed at `date +"%T"`\n" ;
echo
echo -e "ChIP pipeline finished `date +"%T"`"
echo
exit ;
