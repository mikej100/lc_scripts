#!/bin/bash
#SBATCH --partition=long
#SBATCH --ntasks=5
#SBATCH --mem=50G
#SBATCH --time=00-50:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

###########################################################################################################################################################
#	Pipeline for full ChIPmentation analysis													  #
#	Coded by Joseph Blayney                                                                                                                           #
#       Edited by Lucy Cornell Feb 2024														          #
###########################################################################################################################################################

#TO RUN THIS SCRIPT: sbatch ChIPmentation_pipeline_full.sh <genome_name>

alias echo="echo -e"

#Create a date and timestamp for when analysis began
Start_time=`date`
Time=`date +"%T"`
workingdir="$(pwd)"

#Set the genome you are aligning to to $genome, and enter this as the first argument.
genome="$1"

# Designate file name as model name
model=$(echo $(pwd) | awk -F/ '{print $NF}') ;

#Initialise the required modules/packages for analysis.
module unload python-base
module load bowtie2 bedtools/2.25.0 samtools trim_galore flash seqkit ;

#Generally, genome will be mm9 or hg19. Can use if statement to set correct species based on genome build.
#While working on stopgap, can't point towards indices located in t1-data. Therefore have hashed out this if loop, and have
#pointed specifically to the mm9 genome index, built in stopgap. Once we move back to t1-data, will reinstate this chunk, so that
#you can specify which genome to align to as an argument.

if [ "$genome" == "mm9" ] || [ "$genome" == "mm10" ] || [ "$genome" == "mm39" ]; then 
        species="Mus_musculus"
elif [ "$genome" == "hg19" ]; then 
        species="Homo_sapiens"
else
	echo "ChIPmentation pipeline: Please check your genome build and/or species."  
	exit
fi ;


echo "----------------- Bowtie2 Alignment for ChIPmentation----------------- " ;
echo
echo "Analysis started: $Start_time" ;
echo
echo "Current working directory is: "$workingdir" "
echo
echo "Genome Bowtie2 index: "$species" "$genome" "
echo
echo "Sample: "${model}" "
echo

## Print inputs
echo -e "Input fastqs:\n$(ls *.fastq.gz | xargs -n1)" ;
echo

#Concatenate 4XR1 and 4XR2 fastq files into a single R1 and a single R2 file.
nohup cat *L001_R1* *L002_R1* *L003_R1* *L004_R1* > READ1.fastq.gz 
nohup cat *L001_R2* *L002_R2* *L003_R2* *L004_R2* > READ2.fastq.gz 

if [ -s READ1.fastq.gz ]
	then
	echo "READ1 reads successfully concatenated at `date +"%T"`."
        echo
	else
	echo "ChIPmentation pipeline: READ1 concatenation unsuccessfull."  
	exit
fi

if [ -s READ2.fastq.gz ]
        then
        echo "READ2 reads successfully concatenated at `date +"%T"`."
        echo
        else
        echo "ChIPmentation pipeline: READ2 concatenation unsuccessfull."  
        exit
fi

#Use bowtie2 to align the fastq files to the designated genome (set as $genome).
#Currently, allow reads which multi-map to max 2 regions to continue into downstream analysis. Change -k to alter this. 
#Currently, allow maximum of one mismatch in seed. Change -N to alter this. Allowing fewer mismatches will increase accuracy; allowing more 
#will increase sensitivity.
#Currently, maximum fragment length allowed is 1000bp. Change --maxins if you want to alter this e.g. perhaps change to 150bp if you want 
#only mononucleosomes.

nohup bowtie2 -x /databank/igenomes/"$species"/UCSC/"$genome"/Sequence/Bowtie2Index/genome -k 2 -N 1 -1 READ1.fastq.gz -2 READ2.fastq.gz -p 4 --maxins 1000 --un-conc unaligned_1 --al-conc aligned_1 -S alignment.sam 

if [ -s alignment.sam ]
        then
        echo "Alignment successful at `date +"%T"`." 
        echo "See results in alignment.sam. Aligning reads are stored in aligned_1 fastq; unaligning reads stored in unaligned_1 fastq"
        echo
        else
        echo "ChIPmentation pipeline: Alignment unsuccessful."  
        exit
fi

rm -rf *.fastq
rm -rf *.sam

#Trim all reads which did not align in the first bowtie2 alignment.

nohup trim_galore --paired --dont_gzip --length 10 --quality 10 unaligned_1.1 unaligned_1.2

if [ -s unaligned_1.1_val_1.fq ]
        then
        echo "Trimming successful at `date +"%T"`. "
        echo
        else
        echo "ChIPmentation pipeline: Trimming unsuccessful."   
        exit
fi

rm -rf unaligned_1.1 unaligned_1.2

#Retry alignment with the reads which have now been trimmed.

nohup bowtie2 -x /databank/igenomes/"$species"/UCSC/"$genome"/Sequence/Bowtie2Index/genome -k 2 -N 1 -1 unaligned_1.1_val_1.fq -2 unaligned_1.2_val_2.fq -p 4 --maxins 1000 --al-conc aligned_2 --un-conc unaligned_2 -S alignment2.sam 

if [ -s alignment2.sam ]
        then
        echo "Second alignment successful at `date +"%T"`."
        echo
        else
        echo "ChIPmentation pipeline: Second alignment unsuccessful."  
        exit
fi

rm -rf unaligned_1.1_val_1.fq unaligned_1.2_val_2.fq
rm -rf *.sam

#Can now try to flash reads which still don't align. Flashing entails looking for overlap between the mate 1 and mate 2 reads. If there is significant 
#overlap, then it merges the two reads, to give you a new, extended (flashed) read, which is essentially a single end fastq read entry.

nohup flash unaligned_2.1 unaligned_2.2 

if [ -s out.extendedFrags.fastq ]
        then
        echo "Flashing successful at `date +"%T"`."
        echo
        else
        echo "ChIPmentation pipeline: Flashing unsuccessful."   
        exit
fi

rm -rf unaligned_2.1 unaligned_2.2
rm -rf *.sam

#Before flashed reads can be aligned to the genome, need to create a second mate (we want to align paired end reads, so that the resultant files
#are compatible with the alignments conducted up to now). To do this, we generate a fastq file with the reverse complement of the flashed reads file.

nohup seqkit seq -r -p -t dna out.extendedFrags.fastq > out.extendedFragsRevComp.fastq

if [ -s out.extendedFragsRevComp.fastq ]
        then
        echo "Reverse complementation of flashed reads successful at `date +"%T"`."
        echo
        else
        echo "ChIPmentation pipeline: Reverse complementation unsuccessful."   
        exit
fi

#The flashed reads and reverse complement files have exactly the same start and end coordinates. Bowtie2 doesn't like this and won't align them as we want.
#To get around this, we can chew a single bp from the 3' end of each read in each file, to generate 1bp 5' overhangs, allowing bowtie2 to recognise the
#two reads as mate pairs.

nohup cutadapt -u -1 out.extendedFrags.fastq > Flashed_R1.fastq
nohup cutadapt -u -1 out.extendedFragsRevComp.fastq > Flashed_R2.fastq

if [ -s Flashed_R1.fastq ]
        then
        echo "Generation of 1bp overhang successful at `date +"%T"`."
        echo
        else
        echo "ChIPmentation pipeline: Generation of 1bp overhang unsuccessful."   
        exit
fi

rm -rf out.extendedFrags.fastq out.extendedFragsRevComp.fastq
rm -rf *.sam

nohup bowtie2 -x /databank/igenomes/"$species"/UCSC/"$genome"/Sequence/Bowtie2Index/genome -k 2 -N 1 -1 Flashed_R1.fastq -2 Flashed_R2.fastq -p 4 --maxins 1000 --al-conc aligned_3 --un-conc unaligned_3 -S alignment3.sam

if [ -s alignment3.sam ]
        then
        echo "Alignment of flashed reads successful at `date +"%T"`."
        echo
        else
        echo "ChIPmentation pipeline: Alignment of flashed reads unsuccessful."   
        exit
fi

rm -rf Flashed_R1.fastq Flashed_R2.fastq
rm -rf *.sam

#Concatenate all successfully aligned reads into individual READ1 and READ2 fastqs.

cat aligned_*.1 > READ1_al.fastq
cat aligned_*.2 > READ2_al.fastq

if [ -s READ1_al.fastq ]
        then
        echo "Aligning read fastq concatenation successful at `date +"%T"`."
        echo
        else
        echo "ChIPmentation pipeline: Aligning read fastq concatenation unsuccessful."   
        exit
fi

rm -rf *_1.1 *_1.2 *_2.1 *_2.2 *_3.1 *_3.2
rm -rf *.sam

#Do a fourth and final bowtie2 alignment with all aligning reads as input. Alignment rate should be very close to 100%!

nohup bowtie2 -x /databank/igenomes/"$species"/UCSC/"$genome"/Sequence/Bowtie2Index/genome -k 2 -N 1 -1 READ1_al.fastq -2 READ2_al.fastq -p 4 --maxins 1000 -S final_alignment.sam 

if [ -s final_alignment.sam ]
        then
        echo "Final alignment successful at `date +"%T"`."
        echo "See results in final_alignment.sam"
        echo
        else
        echo "ChIPmentation pipeline: Final alignment unsuccessful."  
        exit
fi

#remove intermediate, unzipped fastq files.

rm -rf *.fastq *.fq 

echo -e "Fastq files removed" 
echo

#Filter the mapped reads (filtering the mapped reads from unmapped.

nohup samtools view -@ 8 -bS -F 4 -o mapped.bam final_alignment.sam 

if [ -s mapped.bam ]
        then
        echo "Read filtering successful at `date +"%T"`. "
        echo "See results in mapped.bam."
        echo
        else
        echo "ChIPmentation pipeline: Read filtering unsuccessful."   
        exit
fi

rm -rf *.sam

#Sort the reads in mapped.bam by coordinate.

nohup samtools sort -@ 8 mapped.bam -o mapped.srtC 

if [ -s mapped.srtC ]
        then
        echo "Read sorting successful at `date +"%T"`. "
        echo "See results in mapped.srtC."
        echo
        else
        echo "ChIPmentation pipeline: Read sorting unsuccessful."  
        exit
fi

rm -rf mapped.bam

#Remove PCR duplicates.

nohup samtools rmdup mapped.srtC filtered.bam 

if [ -s filtered.bam ]
        then
        echo "PCR duplicates removal successful at `date +"%T"`. "
        echo "See results in filtered.bam."
        echo
        else
        echo "ChIPmentation pipeline: PCR duplicates removal unsuccessful."   
        exit
fi

rm -rf mapped.srtC *.txt *.histogram *.hist

#Index the bam file

nohup samtools index filtered.bam 

if [ -s filtered.bam.bai ]
        then
        echo "filtered.bam successfully indexed at `date +"%T"`. "
        echo "See results in filtered.bam.bai."
        echo
        else
        echo "ChIPmentation pipeline: bam indexing unsuccessful."   
        exit
fi

module unload python-base
module load python-cbrg

#Create the bigwig (CPM (can swap CPM for RPKM here if you want. Other normalisation options available too) and smoothed over 300bp in 
#50bp bins. 50bp bins is default but can change by altering -bs flag. Can change smoothing window by changing --smoothLength. Can also change 
#normalisation method by altering --normalizeUsing. Increase -bs and --smoothLength metrics for smoother data (less noisy),
#Decrease both to improve resolution...better if sequencing is deeper.

nohup bamCoverage -b filtered.bam --normalizeUsing RPKM -bs 1 -p 4 -o filtered.rpkm.bw 

if [ -s filtered.rpkm.bw ]
        then
        echo "RPKM normalized bigwig successfully generated at `date +"%T"`." 
        echo "See results in filtered.rpkm.bw."
        echo
        else
        echo "ChIPmentation pipeline: bigwig generation unsuccessful."  
        exit
fi


#Get rid of intermediate files to save space. But KEEP filtered.bam for Macs2 analysis (need bam for peak calling and downstream analysis).

rm -rf  mapped.srtC mapped.bam
rm alignment.sam ;

# Change names
mv filtered.bam ${model}.bam
mv filtered.bam.bai ${model}.bam.bai
mv filtered.rpkm.bw ${model}.rpkm.bw

echo -e "Intermediate files removed at `date +"%T"`." 
echo
echo -e "ChIPmentation pipeline finished `date +"%T"`" 
exit 
