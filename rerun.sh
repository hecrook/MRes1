#!/bin/sh
#Hannah Crook
#small script for chnaging names of second runs

#PBS -l select=1:ncpus=4:mem=20gb
#PBS -l walltime=12:0:0
#PBS -N changingnames
#PBS -J 1-26

cd /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/reads2

sample=`ls | tail -n $PBS_ARRAY_INDEX | head -1`
echo "Working with sample "$sample
cd $sample

mv "$sample"_merged_trimmed_R1.fastq 2_"$sample"_merged_trimmed_R1.fastq