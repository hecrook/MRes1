#!/bin/sh
# Hannah Crook
# 16/03/2023

# Project: BrainMets
# Data: sWGS
# Notes: Bwamem on forwards reads only (single-end)

#Array job!

#PBS -l select=1:ncpus=16:mem=62gb
#PBS -l walltime=48:0:0
#PBS -N BrainMetssWGSBWAMEM
#PBS -J 1-26

module load anaconda3/personal
source activate oct_sWGS

fasta=/rds/general/user/hec22/projects/mcneish_team_data/live/reference_data/Ensembl_grch37_primary_sm/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa

tmp=/rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/temp
cd /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/reads2

sample=`ls | tail -n $PBS_ARRAY_INDEX | head -1`
echo "Working with sample "$sample
cd $sample

r1="$sample"_merged_trimmed_R1.fastq
outfile1="$sample".sam
echo "Aligning using bwa-mem.."
bwa mem -M -t 16 $fasta $r1 > $outfile1

echo "Converting sam to bam.."
outfile2="$sample".bam
samtools view -S -b $outfile1 > $outfile2

#Remove sam (storage issues)
rm $outfile1

echo "Sorting bam.."
file3="$sample"_sorted.bam
samtools sort $outfile2 -o $file3

#Remove unsorted bam (storage issues)
rm $outfile2

echo "Indexing the sorted bam.."
samtools index $file3

## QUALIMAP
mkdir -p "$sample"_qualimap_results
echo "Running qualimap.."
conda deactivate 
source activate qualimap
JAVA_OPTS='-Xmx55g'
qualimap bamqc -bam $file3 -outdir "$sample"_qualimap_results --paint-chromosome-limits --genome-gc-distr HUMAN --collect-overlap-pairs -outformat HTML --java-mem-size=55G
