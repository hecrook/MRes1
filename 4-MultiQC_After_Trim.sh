#PBS -lselect=1:ncpus=1:mem=20gb 
#PBS -lwalltime=24:00:00 

module load anaconda3/personal
source activate oct_sWGS 

cd /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/reads

# multiqc . -o /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/MultiQC --ignore *001* --ignore *paired* --filename MultiQC_forwardonly_rerun
multiqc . -o /rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/MultiQC --ignore *R1* --ignore *trimmed* --ignore *merged* --ignore *sorted* --ignore *qualimap* --filename MultiQC_R2_beforetrim