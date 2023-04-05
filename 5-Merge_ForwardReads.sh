#!/bin/sh
# Hannah Crook
# 10/03/2023

# project: BrainMets
# data: sWGS
# notes: Merge all forward, trimmed reads from each sample into one file to take forward into downstream analysis

#Array Job!

#PBS -l select=1:ncpus=4:mem=20gb
#PBS -l walltime=12:0:0
#PBS -N BrainMetssWGS_R1_Merge
#PBS -J 1-26

module load anaconda3/personal
source activate oct_sWGS

tmp=/rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/temp
cd /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/reads2

sample=`ls | tail -n $PBS_ARRAY_INDEX | head -1`
echo "Working with sample "$sample
cd $sample

#remove previously made fastqc files
#cd "$sample"_fastQC/
#rm *merged*
rm "$sample"_merged_trimmed_R1.fastq
#cd ..

#copy over the second run
rsync -avP /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/reads/"$sample"/"$sample"_merged_trimmed_R1.fastq /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/reads2/"$sample"/

echo "Merging forward(R1) reads"
cat "$sample"_merged_trimmed_R1.fastq 2_"$sample"_merged_trimmed_R1.fastq > "$sample"_merged2_trimmed_R1.fastq

rm "$sample"_merged_trimmed_R1.fastq
mv "$sample"_merged2_trimmed_R1.fastq "$sample"_merged_trimmed_R1.fastq

#run fastqc on merged files
r1="$sample"_merged_trimmed_R1.fastq
fastqc $r1 -d $tmp

mv *_fastqc* "$sample"_fastQC/