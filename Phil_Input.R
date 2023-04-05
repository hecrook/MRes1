library(tidyverse)
#TEST RUN SEGMENT FILES 
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/ACE/testrun/")
# Adding a column to each segmentfile which gives the sampleID

segments.tsv <- list.files(getwd(),pattern = "_segments.tsv")
sampleID <- gsub("_segments.tsv","",files)
df.ls <- data.frame("Chromosome" = NA, "Start" = NA, "End" = NA, "Segment_Mean2" = NA, "sample"= NA)
df.ls <- df.ls[-1,]

for(i in 3:length(segments.tsv)) {
  df <- read.table(segments.tsv[i], header = T, sep = "\t", stringsAsFactors = F)
  sampleID <- gsub("_segments.tsv","",segments.tsv[i])
  df <- df %>%
    add_column(sample = sampleID)
  df.ls <- df.ls %>%
    add_row(df[,c(1:3,6,10)])
}

colnames(df.ls)[4] <- "segVal"
colnames(df.ls)[2] <- "start"
colnames(df.ls)[3] <- "end"

colnames(df.ls) <- c('chromosome','start','end','segVal','sample')

###############

df <- read.table(segments.tsv[2], header = T, sep = "\t", stringsAsFactors = F)
sampleID <- gsub("_segments.tsv","",segments.tsv[2])
df <- df %>%
  add_column(sample = sampleID)
df.ls <- df.ls %>%
  add_row(df[,c(1:3,6,10)])
