#!/bin/bash
#SBATCH --account=qchu
#SBATCH --job-name='Chip_align_Luly'
#SBATCH --cpus-per-task=30
#SBATCH --time=24:00:00 
#SBATCH --output=Chip_align_Luly_peak_calling.out


module load samtools
module load bedtools
module load HOMER


#peak calling
list=(
m1643t_noDups_down 
m1644t_noDups_down 
m1649t_noDups_down
m1650t_noDups_down 
m1661t_noDups_down 
m1662t_noDups_down 
m1666t_noDups_down 
m1667t_noDups_down 
m1672t_noDups_down 
m1673t_noDups_down 
m1684t_noDups_down 
m1685t_noDups_down 
)
for an  in ${list[@]}; 
do
makeTagDirectory tag${an} ${an}.bam
findPeaks tag${an} -style factor -o Homer_${an}_spike_in_peak_narrow
done

#FOXA1_ChIP all Peak merging
mergePeaks \
Homer_m1643t_noDups_down_spike_in_peak_narrow \
Homer_m1649t_noDups_down_spike_in_peak_narrow \
Homer_m1684t_noDups_down_spike_in_peak_narrow > FOXA1_DKO_spike_in_xout_all

mergePeaks \
Homer_m1661t_noDups_down_spike_in_peak_narrow \
Homer_m1666t_noDups_down_spike_in_peak_narrow \
Homer_m1672t_noDups_down_spike_in_peak_narrow > FOXA1_SKO_spike_in_xout_all

#FOXA1 PP_vs_PF peak overlapping
mergePeaks FOXA1_DKO_spike_in_xout_all FOXA1_SKO_spike_in_xout_all -prefix over

mergePeaks over_FOXA1_DKO_spike_in_xout_all over_FOXA1_DKO_spike_in_xout_all_FOXA1_SKO_spike_in_xout_all > SKO_only_FOXA1_venn.txt
mergePeaks over_FOXA1_SKO_spike_in_xout_all over_FOXA1_DKO_spike_in_xout_all_FOXA1_SKO_spike_in_xout_all > DKO_only_FOXA1_venn.txt

#peak_annotation
pos2bed.pl over_FOXA1_DKO_spike_in_xout_all > DKO_only_FOAX1_spike_in_xout_all.bed
pos2bed.pl over_FOXA1_DKO_spike_in_xout_all_FOXA1_SKO_spike_in_xout_all > share_only_FOXA1_spike_in_xout_all.bed
pos2bed.pl over_FOXA1_SKO_spike_in_xout_all > SKO_only_FOXA1_spike_in_xout_all.bed
annotatePeaks.pl DKO_only_FOAX1_spike_in_xout_all.bed mm10 > anno_DKO_only_FOAX1_spike_in_xout_all.txt
annotatePeaks.pl share_only_FOXA1_spike_in_xout_all.bed mm10 > anno_share_only_FOXA1_spike_in_xout_all.txt
annotatePeaks.pl SKO_only_FOXA1_spike_in_xout_all.bed mm10 > anno_SKO_only_FOXA1_spike_in_xout_all.txt