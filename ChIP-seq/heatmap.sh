#!/bin/bash
#SBATCH --account=qchu
#SBATCH --job-name='Chip_align_Luly'
#SBATCH --cpus-per-task=40
#SBATCH --time=24:00:00 
#SBATCH --output=Chip_align_Luly_heatmap.out

module load deeptools


#Figure 3 C
computeMatrix reference-point \
--referencePoint center \
-a 5000 -b 5000 \
-S /data1/qchu/ChIP_seq/Luly/bw/m1661t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1666t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1672t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1643t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1649t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1684t_noDups_down.bw \
--samplesLabel m1661 m1666 m1672 m1643 m1649 m1684 \
-R SKO_only_FOXA1_spike_in_xout_all.bed \
--missingDataAsZero \
-p 40 \
-o FOXA1_spike_in_over_xinput_xout_all_P_only.matrix \
-bs 20

plotHeatmap -m FOXA1_spike_in_over_xinput_xout_all_P_only.matrix \
--sortRegions descend \
--sortUsing mean --sortUsingSamples 1 2 3 \
--averageTypeSummaryPlot mean \
--startLabel center \
--regionsLabel "SKO_only" \
--zMax 2.5 \
--dpi 400 --colorList 'white,red' \
-o FOXA1_spike_in_over_xinput_xout_all_spike_in.pdf \



#Figure 3 F
mergePeaks -d 10 w18_SKO_only_FOXA1_0.1unchange_peak_highest.bed xFOXA1_0.1_unchange_position_TSS.bed > all_FOXA1_0.1_unchange_position
mergePeaks -d 10 w18_FOXA1_1.0up_peak_highest.bed xFOXA1_1.0_DEG_up_position_TSS.bed > all_FOXA1_1.0_DEG_up_position
mergePeaks -d 10 w18_FOXA1_1.0down_peak_highest.bed xFOXA1_1.0_DEG_down_position_TSS.bed > all_FOXA1_1.0_DEG_down_position

pos2bed.pl all_FOXA1_0.1_unchange_position > all_FOXA1_0.1_unchange_position.bed
pos2bed.pl all_FOXA1_1.0_DEG_up_position > all_FOXA1_1.0_DEG_up_position.bed
pos2bed.pl all_FOXA1_1.0_DEG_down_position > all_FOXA1_1.0_DEG_down_position.bed

computeMatrix reference-point \
--referencePoint center \
-a 5000 -b 5000 \
-S /data1/qchu/ChIP_seq/Luly/bw/m1661t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1666t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1672t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1643t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1649t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1684t_noDups_down.bw \
--samplesLabel m1661 m1666 m1672 m1643 m1649 m1684 \
-R  all_FOXA1_1.0_DEG_up_position.bed all_FOXA1_0.1_unchange_position.bed all_FOXA1_1.0_DEG_down_position.bed \
--missingDataAsZero \
-p 40 \
--outFileSortedRegions FOXA1_sorted.bed \
-o FOXA1_DEG.matrix \
-bs 20

plotHeatmap -m FOXA1_DEG.matrix \
--sortRegions descend \
--sortUsing mean --sortUsingSamples 1 2 3 \
--averageTypeSummaryPlot mean \
--startLabel center \
--regionsLabel "up" "unchange" "down" \
--zMax 3 --zMin 0 \
--dpi 400 --colorList 'white,red' \
--whatToShow 'heatmap and colorbar' \
-o FOXA1_DEG_heatmap.pdf \


plotProfile -m FOXA1_DEG.matrix \
--averageType mean \
--startLabel center \
--regionsLabel "up" "unchange" "down" \
-o FOXA1_DEG_density.pdf \
--dpi 400 --colors '#A41125' '#8C9599' '#3482B6'



#Figure 3 G
computeMatrix reference-point \
--referencePoint center \
-a 5000 -b 5000 \
-S /data1/qchu/ChIP_seq/Luly/bw/m1662t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1667t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1673t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1644t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1650t_noDups_down.bw \
/data1/qchu/ChIP_seq/Luly/bw/m1685t_noDups_down.bw \
--samplesLabel m1662 m1667 m1673 m1644 m1650 m1685 \
-R FOXA1_order_up_TSS_FOXA1bind.bed FOXA1_order_unchange_TSS_FOXA1bind.bed FOXA1_order_down_TSS_FOXA1bind.bed \
--missingDataAsZero \
-p 40 \
-o k27ac_on_FOXA1_gene_tss_FOXA1bind.matrix \
-bs 20

plotHeatmap -m k27ac_on_FOXA1_gene_tss_FOXA1bind.matrix \
--sortRegions keep \
--averageTypeSummaryPlot mean \
--startLabel center \
--regionsLabel "up" "unchange" "down" \
--zMax 30 --zMin 0 \
--dpi 400 --colorList 'white,red' \
--whatToShow 'heatmap and colorbar' \
-o k27ac_on_FOXA1_gene_tss_FOXA1bind.pdf \


plotProfile -m k27ac_on_FOXA1_gene_tss_FOXA1bind.matrix \
--averageType mean \
--startLabel center \
--regionsLabel "up" "unchange" "down" \
-o FOXA1_DEG_density.pdf \
--dpi 400 --colors '#A41125' '#8C9599' '#3482B6'



