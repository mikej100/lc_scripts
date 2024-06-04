#!/bin/bash
#SBATCH --partition=long
#SBATCH --ntasks=5
#SBATCH --mem=50G
#SBATCH --time=00-50:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err
#SBATCH --mail-user=lucy.cornell@imm.ox.ac.uk
#SBATCH --mail-type=end,fail

module load python-cbrg bedtools

cd /project/higgslab/lcornell/ChIPs/8_CD71_EB_Rad21_INVEN/Peak_analysis/Lanceotron_outputs/

for bed in CD71*_toppeaks.bed ; do 
    echo "Sorting $bed" ;
    bedtools sort -i "$bed" > sorted_"$bed" ;
done

# Intersecting postions from E14 WT and INVEN samples (A4.2) - B12 was removed as it has a HS38 deletion
echo "Finding shared peaks across the samples" ;

if ! bedtools intersect -a sorted_CD71_EB_Rad21_Mira_INVEN_E14_Rep1_L-tron_toppeaks.bed \
 -b sorted_CD71_EB_Rad21_Mira_INVEN_E14_Rep2_L-tron_toppeaks.bed \
 sorted_CD71_EB_Rad21_Mira_INVEN_A42_Rep1_L-tron_toppeaks.bed \
 sorted_CD71_EB_Rad21_Mira_INVEN_A42_Rep2_L-tron_toppeaks.bed \
 -sorted > temp_Rad21_shared_toppeaks.bed ; then
    echo "Bedtools intersect failed" ;
    exit 1 ;
fi

echo "Merging shared peaks so there aren't any repeated regions" ;

awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$8}' temp_Rad21_shared_toppeaks.bed | sort -k1,1 -k2,2n > sorted_temp_Rad21_shared_toppeaks.bed

if ! bedtools merge -i sorted_temp_Rad21_shared_toppeaks.bed -c 4 -o mean | awk 'BEGIN{FS=OFS="|"} {print $1}'> Rad21_shared_toppeaks.bed ; then
    echo "Bedtools merge failed" ;
    exit 1 ;
fi

echo "Starting multicov" ;

cd /project/higgslab/lcornell/ChIPs/8_CD71_EB_Rad21_INVEN/

if ! bedtools multicov -bams CD71_EB_Rad21_Mira_INVEN_E14_Rep1/CD71_EB_Rad21_Mira_INVEN_E14_Rep1.bam \
 CD71_EB_Rad21_Mira_INVEN_E14_Rep2/CD71_EB_Rad21_Mira_INVEN_E14_Rep2.bam \
 CD71_EB_Rad21_Mira_INVEN_A4.2_Rep1/CD71_EB_Rad21_Mira_INVEN_A4.2_Rep1.bam \
 CD71_EB_Rad21_Mira_INVEN_A4.2_Rep2/CD71_EB_Rad21_Mira_INVEN_A4.2_Rep2.bam \
 CD71_EB_Rad21_Mira_INVEN_B12_Rep1/CD71_EB_Rad21_Mira_INVEN_B12_Rep1.bam \
 CD71_EB_Rad21_Mira_INVEN_B12_Rep2/CD71_EB_Rad21_Mira_INVEN_B12_Rep2.bam \
 -bed Peak_analysis/Lanceotron_outputs/Rad21_shared_toppeaks.bed > Peak_analysis/Lanceotron_outputs/temp_Coverage_over_shared_toppeaks.bed ; then 
    echo "Multicov failed" ;
    exit 1;
fi

echo -e "chrom\tstart\tend\tregion_id\tE14_r1\tE14_r2\tA42_r1\tA42_r2\tB12_r1\tB12_r2" > header ;

cat header Peak_analysis/Lanceotron_outputs/temp_Coverage_over_shared_toppeaks.bed > Peak_analysis/Lanceotron_outputs/Coverage_over_shared_toppeaks.bed ;

rm header Peak_analysis/Lanceotron_outputs/*temp* ;

# Pulling peaks from the results for comparisons
cd /project/higgslab/lcornell/ChIPs/8_CD71_EB_Rad21_INVEN/Peak_analysis/Lanceotron_outputs/

##### Hba 
# 34690 = Snrp23 promoter
# 34691 = Rhbdf1 CTCF
# 37570 = Rhbdf1 promoter
# 34693 = Mpg promoter
# 37571 = HS38
# 37572 = HS39
# 37574 = R2 enhancer
# 37582 = Hba gene 

##### Hbb
# 24298 = CTCF
# 24301 = CTCF_2
# 25698 = enhancer
# 25701 = CTCF_enhancer

echo "Complete!" ;

