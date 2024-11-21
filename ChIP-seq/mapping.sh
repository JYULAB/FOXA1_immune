#!/bin/bash
#SBATCH --account=qchu
#SBATCH --job-name='Chip_align_Luly'
#SBATCH --cpus-per-task=40
#SBATCH --time=24:00:00 
#SBATCH --output=FOXA1_H3k27ac_Chip_align.out

module load bowtie2
module load samtools
module load picard
module load Trimmomatic

list=(
m1643
m1644
m1649
m1650
m1661
m1662
m1666
m1667
m1672
m1673
m1684
m1685
)
for an  in ${list[@]}; 
do
#mm10 mapping
java -jar /home/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 40 ${an}_R1.fastq ${an}_R2.fastq ${an}t_R1.fastq ${an}t_R1_unpaired_trimmed.fastq ${an}t_R2.fastq ${an}t_R2_unpaired_trimmed.fastq ILLUMINACLIP:/home/apps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
bowtie2 -x /home/apps/iGenomes/references/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 ${an}t_R1.fastq -2 ${an}t_R2.fastq -S ${an}t.sam
samtools view -b -S ${an}t.sam > ${an}t.bam
samtools sort -o ${an}t_sorted.bam ${an}t.bam
samtools index ${an}t_sorted.bam
java -Xmx1G -jar /home/apps/picard-3.0.0/picard.jar MarkDuplicates I=${an}t_sorted.bam O=${an}t_noDups.bam M=${an}t_noDups.txt REMOVE_DUPLICATES=true

#Drosophila spike_in mapping
bowtie2 -x /data1/czhao/xiaodong/p300/Drosophila/UCSC/dm6/Sequence/Bowtie2Index/genome -1 ${an}t_R1.fastq -2 ${an}t_R2.fastq -S ${an}t_drosoph.sam
samtools view -b -S ${an}t_drosoph.sam > ${an}t_drosoph.bam
samtools sort -o ${an}t_drosoph_sorted.bam ${an}t_drosoph.bam
samtools index ${an}t_drosoph_sorted.bam
java -Xmx1G -jar /home/apps/picard-3.0.0/picard.jar MarkDuplicates I=${an}t_drosoph_sorted.bam O=${an}t_drosoph_noDups.bam M=${an}t_drosoph_noDups.txt REMOVE_DUPLICATES=true
done