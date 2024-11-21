#!/bin/bash
#SBATCH --account=qchu
#SBATCH --job-name='Chip_align_Luly'
#SBATCH --cpus-per-task=50
#SBATCH --time=48:00:00 
#SBATCH --output=Luly_spike_in_down_sampling.out

module load bowtie2
module load samtools
module load bedtools
module load deeptools

#FOXA1_ChIP_down_sampling
samtools view -@ 40 -s 0.8624 -b m1643t_noDups.bam -o m1643t_noDups_down.bam


samtools view -@ 40 -s 1.00 -b m1649t_noDups.bam -o m1649t_noDups_down.bam


samtools view -@ 40 -s 0.2115 -b m1661t_noDups.bam -o m1661t_noDups_down.bam


samtools view -@ 40 -s 0.0940 -b m1666t_noDups.bam -o m1666t_noDups_down.bam


samtools view -@ 40 -s 0.1121 -b m1672t_noDups.bam -o m1672t_noDups_down.bam


samtools view -@ 40 -s 0.2440 -b m1684t_noDups.bam -o m1684t_noDups_down.bam


#H3K27ac_ChIP_down_sampling
samtools view -@ 40 -s 1.00 -b m1644t_noDups.bam -o m1644t_noDups_down.bam


samtools view -@ 40 -s 0.7078 -b m1650t_noDups.bam -o m1650t_noDups_down.bam


samtools view -@ 40 -s 0.6105 -b m1662t_noDups.bam -o m1662t_noDups_down.bam


samtools view -@ 40 -s 0.2244 -b m1667t_noDups.bam -o m1667t_noDups_down.bam


samtools view -@ 40 -s 0.2276 -b m1673t_noDups.bam -o m1673t_noDups_down.bam


samtools view -@ 40 -s 0.2408 -b m1685t_noDups.bam -o m1685t_noDups_down.bam


#bigWig file generate:
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
samtools index -@40 ${an}.bam
bamCoverage --bam ${an}.bam -o ${an}.bw \
--binSize 10 --smoothLength 10 --normalizeUsing None --numberOfProcessors 40
done




