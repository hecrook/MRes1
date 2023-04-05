#!/bin/sh
# Hannah Crook
# 10/03/2023

# Project: BrainMets_sWGS
# Data type: sWGS
# Notes: Run trimmomatic on raw sWGS fastq files

#Array Job!

#PBS -l select=1:ncpus=4:mem=20gb
#PBS -l walltime=12:0:0
#PBS -N OCTOPUS-Trim
#PBS -J 1-26

module load anaconda3/personal
source activate oct_sWGS
# source activate picard

tmp=/rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/temp
af=/rds/general/user/hec22/projects/mcneish_team_data/live/Hannah_CNS/analyses/Script/adapters.fa

cd /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/reads2

sample=`ls | tail -n $PBS_ARRAY_INDEX | head -1`
echo "Working with sample "$sample
cd $sample

#Remove merged and merged_trimmed fastQC files
#cd "$sample"_fastQC
#rm *merged_R*_fastqc.*
#rm *trimmed*
#rm *trimmed*

# Remove merged and trimmed files
#cd ..
#rm *merged*
#rm *trimmed*

# Moving FastQC results into a separate folder (declutter)
mkdir -p "$sample"_fastQC
mv *_fastqc* "$sample"_fastQC/

# Unzipping all files (required for merging lanes)
gunzip -r *.gz .


# Merging fastq files across lanes (This cannot be done at this point IF THERE ARE MULTIPLE RUNS as we need to add read groups before the files are merged. Addreadgroups must be done on a sam or bam file, meaning we must have already trimmed and aligned our fastq files. Therefore we must do trimmomatic, and bwamem on each individual run of the same sample before we can add read groups, and then we can merge the files after this point)
#echo "Merging forward(R1) reads"
#cat $(ls | grep R1_001.fastq$) > "$sample"_merged_R1.fastq
#echo "Merging reverse(R2) reads"
#cat $(ls | grep R2_001.fastq$) > "$sample"_merged_R2.fastq

# Running FastQC on merged files
#find . -name "*merged_R1.fastq" | xargs -n 1 fastqc -d $tmp
#find . -name "*merged_R2.fastq" | xargs -n 1 fastqc -d $tmp


#r1="$sample"_merged_R1.fastq
#r2="$sample"_merged_R2.fastq

#rename fastq files so they all have the same format
mv *R1_001.fastq "$sample"_R1.fastq
mv *R2_001.fastq "$sample"_R2.fastq

r1="$sample"_R1.fastq
r2="$sample"_R2.fastq

outfile1="$sample"_R1_trimmed_paired.fastq
outfile2="$sample"_R1_trimmed_unpaired.fastq
outfile3="$sample"_R2_trimmed_paired.fastq
outfile4="$sample"_R2_trimmed_unpaired.fastq

# Trimming
trimmomatic PE -threads 4 -phred33 \
			$r1 $r2 \
			$outfile1 $outfile2 \
			$outfile3 $outfile4 \
			ILLUMINACLIP:$af:2:30:10 \
			MINLEN:50

# FastQC on output files
fastqc $outfile1 -d $tmp
fastqc $outfile2 -d $tmp
fastqc $outfile3 -d $tmp
fastqc $outfile4 -d $tmp

mv *_fastqc* "$sample"_fastQC/
