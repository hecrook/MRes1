#!/bin/sh
# Hannah Crook
# 16/03/2023
# Project: BrainMets_sWGS_rerun
# Data type: sWGS
# Notes: Run MultiQC on FastQC files

#PBS -lselect=1:ncpus=1:mem=20gb 
#PBS -lwalltime=24:00:00 
#PBS -N MultiQC-OHA-AmpliSeq

module load anaconda3/personal
source activate oct_sWGS

cd /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/reads2
multiqc . -o /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/MultiQC --filename MultiQC_2ndRun