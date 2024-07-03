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

###############################################################################
#                         Get command line options
# Function to display help text
usage() {
    echo "Utility to extend regions in .bed file by given number of base pairs"
    echo "Usage: $(basename $0) -s <bp> [-gV]"
    echo "Options:"
    echo "  -g          Genome alignment. Default is mm39"
    echo "  -h          Display this help message"
    echo "  -s          Slop distance: bp to add to each end of region"
    echo "  -V          very verbose: print every command line"
}
# Defaults
endedness="paired"
genome="mm39"

# Parse options using getopts
while getopts "hb:V" option; do
    case "${option}" in
        g)  genome=($OPTARG)
            ;;
        b)  slop=($OPTARG)
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
echo "  slop: ${slop}"
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
# Exit script on any error
set -e

#Model name is the folder name
model="$(basename $PWD)"

if [ "$genome" == "mm9" ] || [ "$genome" == "mm10" ] || [ "$genome" == "mm39" ]; then 
        species="Mus_musculus"
elif [ "$genome" == "hg19" ] || [ "$genome" == "hg38" ] ; then 
        species="Homo_sapiens"
else
	echo "Please check your genome build and/or species."  
	exit
fi
ref_genome="/databank/igenomes/"$species"/UCSC/"$genome"/Sequence/Bowtie2Index/genome"
chr_sizes="/databank/igenomes/"$species"/UCSC/"$genome"/Sequence/WholeGenomeFasta/chr_sizes.txt"
echo "chr_sizes: ${chr_sizes}"
#/databank/igenomes/Mus_musculus/UCSC/mm39/Sequence/WholeGenomeFasta/chr_sizes.txt
################################################################################
#                           Set up modules 
module purge #Unloads all modules from the users environment
module load bedtools2/2.27.1
module load ucsctools/385

################################################################################
#                       Run process steps
#                                                                              #
# Run bedtools slop to increase size of regions for subsequent coverage counts
# Expand in both directions by given number of base pairs.

inputf="Lanceotron_output/${model}_L-tron_toppeaks"

echo "Calling bedtools slop for ${inputf}.bed"
bedtools slop \
        -i ${inputf}.bed \
        -g ${chr_sizes} \
        -b ${slop} \
        > "${inputf}_slop${slop}.bed"

echo "$(now) Starting sort for ${inputf}_slop${slop}.bed"
sort -k1,1 -k2,2n \
	"${inputf}_slop${slop}.bed" > \
	"${inputf}_slop${slop}.sorted.bed"

echo "$(now) Starting bedToBigBed for ${inputf}_slop${slop}.sorted.bed$"
bedToBigBed \
	${inputf}_slop${slop}.sorted.bed \
 	$chr_sizes \
	${inputf}_slop${slop}.bb 
echo
echo -e "$(now) Process steps complete"
echo
exit ;
