#!/bin/bash
#SBATCH --partition=short
#SBATCH --ntasks=5
#SBATCH --time=00-05:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err
#SBATCH --mail-type=end,fail

module load python-cbrg bedtools

## mm39 
# chr11:32148464-32175527 = 5' regions covering Snrnp25 and Rhbdf1
# chr11:32224569-32251632 = 3' region covering Hba-x to Hba-a2
# chr7:103479350-103506413 = 3' region of Hbb covering Hbb-bh2 to Hbb-y

# chr11:32188452–32249902 = Emily's defined region
# chr11:32323049–32384895 = Emily's defined region

# Double checking sizes
# awk 'BEGIN{FS=OFS="\t"} {print $0, ($3-$2)}' Manually_defined_regions.bed 

cd /project/higgslab/lcornell/ChIPs/8_CD71_EB_Rad21_INVEN/

if ! bedtools multicov -bams CD71_EB_Rad21_Mira_INVEN_E14_Rep1/CD71_EB_Rad21_Mira_INVEN_E14_Rep1.bam \
 CD71_EB_Rad21_Mira_INVEN_E14_Rep2/CD71_EB_Rad21_Mira_INVEN_E14_Rep2.bam \
 CD71_EB_Rad21_Mira_INVEN_A4.2_Rep1/CD71_EB_Rad21_Mira_INVEN_A4.2_Rep1.bam \
 CD71_EB_Rad21_Mira_INVEN_A4.2_Rep2/CD71_EB_Rad21_Mira_INVEN_A4.2_Rep2.bam \
 CD71_EB_Rad21_Mira_INVEN_B12_Rep1/CD71_EB_Rad21_Mira_INVEN_B12_Rep1.bam \
 CD71_EB_Rad21_Mira_INVEN_B12_Rep2/CD71_EB_Rad21_Mira_INVEN_B12_Rep2.bam \
 -bed Peak_analysis/Manually_defined_regions.bed > Peak_analysis/temp_Coverage_over_defined_regions.bed ; then 
    echo "Multicov failed" ;
    exit 1;
fi

echo -e "chrom\tstart\tend\tregion_id\tE14_r1\tE14_r2\tA42_r1\tA42_r2\tB12_r1\tB12_r2" > header ;

cat header Peak_analysis/temp_Coverage_over_defined_regions.bed > Peak_analysis/Coverage_over_defined_regions.bed ;

rm header Peak_analysis/*temp* ;

echo "Complete!" ;