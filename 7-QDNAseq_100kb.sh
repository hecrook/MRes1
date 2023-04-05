# Hannah Crook
# 13/03/2023

# Project: BrainMets
# Data: sWGS
# Notes: Running QDNASeq on sWGS forward reads ONLY from BrainMets project

#PBS -lselect=1:ncpus=4:mem=20gb
#PBS -lwalltime=12:0:0

export PATH=/rds/general/user/hec22/home/anaconda3/bin:$PATH
cd /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/scripts
Rscript 7-QDNAseq_100kb.R