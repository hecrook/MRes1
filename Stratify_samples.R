######STRATIFY SAMPLES ON SE Vs ARC SAMPLES#####
#To be done after sigAct has been fetched from CINSignatures
# load data in
meta <- read.table("metadata_HM.txt", header = TRUE, sep = "\t", na.strings = "")
meta <- meta[-which(is.na(meta$PatientNumber)),]
head(meta)
patient <- read.table("Patients_matching_anotherone_HM.txt", header = TRUE, sep = "\t")
head(patient)
####Have the IGF IDs for all of the Ampliseq, but not for the sWGS files, need these to stratify the patients based on study entry and archival. Will do on the subset for now which I have the matched IGFs for.

#Load data: subset IGF data
getwd()
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/Phil_CNS/")
IGF.sWGS <- read.csv("IGF_OCTOPUS_subset.csv", header = TRUE, sep = ",")
head(IGF.sWGS)

#Make 2 different objects, 1 for archival, and 1 for studyentry
arc.b <- IGF.sWGS$SampleType == "Archival"
arc.IGF <- IGF.sWGS[arc.b.2,]

se.b <- IGF.sWGS$SampleType == "StudyEntry"
se.IGF <- IGF.sWGS[se.b,]

sigAct.df <- as.data.frame(sigAct)
sig.arc <- sigAct.df[intersect(rownames(sigAct.df), arc.IGF$IGF_sWGS),]
sig.se <- sigAct.df[intersect(rownames(sigAct.df), se.IGF$IGF_sWGS),]

#####save outputs as .csv files#####
write.csv(sig.arc, file = 'C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/Phil_CNS/sigarc_Drews_OHA', row.names = TRUE)
write.csv(sig.se, file = 'C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/Phil_CNS/sigse_Drews_OHA', row.names = TRUE)
