library(tidyverse)
#TEST RUN SEGMENT FILES 
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/CNS_fig/data/OCTOPUS_segmentfiles/")
# Adding a column to each segmentfile which gives the sampleID

segments.tsv <- list.files(getwd(),pattern = "_segments.tsv")
sampleID <- gsub("_segments.tsv","",segments.tsv)
df.ls <- data.frame("Chromosome" = NA, "Start" = NA, "End" = NA, "Segment_Mean2" = NA, "sample"= NA)
df.ls <- df.ls[-1,]

df.ls$sample <- as.character(df.ls$sample)

for(i in 1:length(segments.tsv)) {
  df <- read.table(segments.tsv[i], header = T, sep = "\t", stringsAsFactors = F)
  sampleID <- gsub("_segments.tsv","",segments.tsv[i])
  df <- df %>%
    add_column(sample = sampleID)
  df.ls <- df.ls %>%
    add_row(df[,c(1:3,6,10)])
}

# colnames(df.ls)[4] <- "segVal"
# colnames(df.ls)[2] <- "start"
# colnames(df.ls)[3] <- "end"

colnames(df.ls) <- c('chromosome','start','end','segVal','sample')

setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/CNS_fig/data/")
write.csv(df.ls, file = "BrainMets_Run2_PhilCNS_input.csv")
##################################CNS CALC
library(CINSignatureQuantification)
#need to subset so only primary samples
df.ls.subset <- df.ls[which(df.ls$sample %in% igf.primary),]

cnobj <- quantifyCNSignatures(object = df.ls.subset,
                              experimentName = "BrainMetsRun2CNS",
                              method = "mac",
                              cores = 1,
                              build = "hg19")
cnobj

#Plotting Functions
#Sample can be specified as a numerical index or sample name
plotSegments(object = cnobj, sample = 1, cn.max = 8)

#Additional arguments can be provided for the heatmap() function
plotSampleByComponent(object = cnobj)

#Custom colours for signatures can be provided with the cols argument, the length of which must match the number of signatures
plotActivities(object = cnobj,type = "threshold")

SxC <- getSampleByComponent(cnobj)
NMF::aheatmap(SxC,fontsize = 7,Rowv=FALSE,Colv=FALSE,legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")

#CNS Calc for primary
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/CNS_fig/data/BrainMets_run2_segmentfiles/")
# Adding a column to each segmentfile which gives the sampleID

segments.tsv <- list.files(getwd(),pattern = "_segments.tsv")
sampleID <- gsub("_segments.tsv","",segments.tsv)
df.ls <- data.frame("Chromosome" = NA, "Start" = NA, "End" = NA, "Segment_Mean2" = NA, "sample"= NA)
df.ls <- df.ls[-1,]

df.ls$sample <- as.character(df.ls$sample)

for(i in 1:length(segments.tsv)) {
  df <- read.table(segments.tsv[i], header = T, sep = "\t", stringsAsFactors = F)
  sampleID <- gsub("_segments.tsv","",segments.tsv[i])
  df <- df %>%
    add_column(sample = sampleID)
  df.ls <- df.ls %>%
    add_row(df[,c(1:3,6,10)])
}

# colnames(df.ls)[4] <- "segVal"
# colnames(df.ls)[2] <- "start"
# colnames(df.ls)[3] <- "end"

colnames(df.ls) <- c('chromosome','start','end','segVal','sample')
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/CNS_fig/data/")
write.csv(df.ls, file = "OCTOPUS_PhilCNS_input.csv")
##################################CNS CALC
library(CINSignatureQuantification)
#need to subset so only archival samples?
df.ls.subset <- df.ls[which(df.ls$sample %in% igf.arc),]

cnobj <- quantifyCNSignatures(object = df.ls.subset,
                              experimentName = "OCTOPUSCNS",
                              method = "mac",
                              cores = 1,
                              build = "hg19")
cnobj

#Plotting Functions
#Sample can be specified as a numerical index or sample name
plotSegments(object = cnobj, sample = 5, cn.max = 8)

#Additional arguments can be provided for the heatmap() function
plotSampleByComponent(object = cnobj)

#Custom colours for signatures can be provided with the cols argument, the length of which must match the number of signatures
plotActivities(object = cnobj,type = "threshold")

##################Run1
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/CNS_fig/data/BrainMets_run1_segmentfiles/")
segments.tsv <- list.files(getwd(),pattern = "_segments.tsv")
sampleID <- gsub("_segments.tsv","",segments.tsv)
df.ls <- data.frame("Chromosome" = NA, "Start" = NA, "End" = NA, "Segment_Mean2" = NA, "sample"= NA)
df.ls <- df.ls[-1,]

df.ls$sample <- as.character(df.ls$sample)

for(i in 1:length(segments.tsv)) {
  df <- read.table(segments.tsv[i], header = T, sep = "\t", stringsAsFactors = F)
  sampleID <- gsub("_segments.tsv","",segments.tsv[i])
  df <- df %>%
    add_column(sample = sampleID)
  df.ls <- df.ls %>%
    add_row(df[,c(1:3,6,10)])
}
colnames(df.ls) <- c('chromosome','start','end','segVal','sample')
write.csv(df.ls, file = "Phil_CNS_input.csv")

df.ls.subset <- df.ls[which(df.ls$sample %in% igf.primary),]

cnobj <- quantifyCNSignatures(object = df.ls.subset,
                              experimentName = "BrainMetsRun1CNS",
                              method = "mac",
                              cores = 1,
                              build = "hg19")
cnobj

#Plotting Functions
#Sample can be specified as a numerical index or sample name
plotSegments(object = cnobj, sample = 1, cn.max = 8)

#Additional arguments can be provided for the heatmap() function
plotSampleByComponent(object = cnobj, Rowv = NULL)

#Custom colours for signatures can be provided with the cols argument, the length of which must match the number of signatures
plotActivities(object = cnobj,type = "threshold")

SxC <- getSampleByComponent(cnobj)
NMF::aheatmap(SxC,fontsize = 7, Rowv = NA, Colv = NA, legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")

#reorder SxC plots
getwd()
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/CNS_fig/data/")
list.files(getwd())
df.ls <- read.csv("OCTOPUS_PhilCNS_input.csv", row.names = 1)
cnobj <- quantifyCNSignatures(object = df.ls,
                              experimentName = "OCTOPUSCNS",
                              method = "mac",
                              cores = 1,
                              build = "hg19")
saveRDS(cnobj, "OCTOPUS.cnobj.Rdata")
SxC <- getSampleByComponent(cnobj)
NMF::aheatmap(SxC,fontsize = 7, Rowv = FALSE, Colv = FALSE, legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")
#creating a specific order for te heatmaps
hclust_rows <- as.dendrogram(hclust(dist((SxC))))
hclust_cols <- as.dendrogram(hclust(dist(t(SxC))))
NMF::aheatmap(SxC,fontsize = 7, Rowv = hclust_rows, Colv = NA, legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")

#will just not order the columns, and then use hclust to order the rows for more order 

setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/CNS_fig/data")
list.files(getwd())
df.ls <- read.csv("BrainMets_Run2_PhilCNS_input.csv", row.names = 1)
cnobj <- quantifyCNSignatures(object = df.ls,
                              experimentName = "Run2CNS",
                              method = "mac",
                              cores = 1,
                              build = "hg19")
saveRDS(cnobj, "Run2.cnobj.Rdata")
SxC <- getSampleByComponent(cnobj)
hclust_rows <- as.dendrogram(hclust(dist((SxC))))
NMF::aheatmap(SxC,fontsize = 14, Rowv = hclust_rows, Colv = NA, legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")

#increase font size
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/CNS_fig/cnobj/")
source("Run1.cnobj.Rdata")
df.ls <- read.csv("Phil_CNS_input.csv", row.names = 1)
cnobj <- quantifyCNSignatures(object = df.ls,
                              experimentName = "Run1CNS",
                              method = "mac",
                              cores = 1,
                              build = "hg19")
save(cnobj, file = "Run1.cnobj.Rdata")
SxC <- getSampleByComponent(cnobj)
hclust_rows <- as.dendrogram(hclust(dist((SxC))))
NMF::aheatmap(SxC,fontsize = 14, Rowv = hclust_rows, Colv = NA, legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")

df.ls <- read.csv("OCTOPUS_PhilCNS_input.csv", row.names = 1)
cnobj <- quantifyCNSignatures(object = df.ls,
                              experimentName = "OCTOPUSCNS",
                              method = "mac",
                              cores = 1,
                              build = "hg19")
save(cnobj, file = "OCTOPUS.cnobj.Rdata")
SxC <- getSampleByComponent(cnobj)
hclust_rows <- as.dendrogram(hclust(dist((SxC))))
NMF::aheatmap(SxC,fontsize = 14, Rowv = hclust_rows, Colv = NA, legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")
