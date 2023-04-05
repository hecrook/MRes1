######### CREATE BOXPLOT OF SEQUENCES BEFORE AND AFTER TRIMMING ########
# Need to pull out information from MultiQC general statistics
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/Sequences_fig/MultiQC_R1_beforetrim_data/")
list.files(getwd())
seq.r1.before <- read.delim("multiqc_general_stats.txt", header = TRUE)
seq.r1.before %>% as_tibble()
colnames(seq.r1.before)
seq.r1.before <- seq.r1.before[,c(1,7)]

setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/Sequences_fig/MultiQC_R1_aftertrim_data/")
seq.r1.after <- read.delim("multiqc_general_stats.txt", header = TRUE)
seq.r1.after %>% as_tibble()
seq.r1.after <- seq.r1.after[,c(1,7)]

setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/Sequences_fig/MultiQC_R2_beforetrim_data/")
seq.r2.before <- read.delim("multiqc_general_stats.txt", header = TRUE)
seq.r2.before %>% as_tibble()
seq.r2.before <- seq.r2.before[,c(1,7)]

setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/Sequences_fig/MultiQC_R2_aftertrim_data/")
seq.r2.after <- read.delim("multiqc_general_stats.txt", header = TRUE)
seq.r2.after %>% as_tibble()
seq.r2.after <- seq.r2.after[,c(1,7)]

#Create column with total sequences for both paired and unpaired
# seq.r1.before has not been trimmed therefore does not have paired and unapired sequences, only total, so needs to be managed differently
seq.r1.before[1] <- gsub("(IGF[0-9]+)_.*", "\\1", seq.r1.before$Sample)
seq.r1.before <- seq.r1.before %>% 
  add_column(stat)
#seq.r1.after
seq.r1.after %>% as_tibble()
sampleID = gsub("(IGF[0-9]+)_.*", "\\1", seq.r1.after$Sample)
seq.r1.after <- seq.r1.after %>%
  add_column(sampleID = sampleID, .before = "Sample")

for (i in 1:length(seq.r1.after$Sample)){
  if (str_detect(seq.r1.after$Sample[i], "unpaired") == TRUE){
    seq.r1.after$Sample[i] <- "unpaired"
  }
  else
    seq.r1.after$Sample[i] <- "paired"
  }


pivot_wider(seq.r1.after, names_from = Sample,
              values_from = FastQC_mqc.generalstats.fastqc.total_sequences)

#seq.r2.before
seq.r2.before %>% as_tibble()
sampleID = gsub("(IGF[0-9]+)_.*", "\\1", seq.r2.before$Sample)
seq.r2.before <- seq.r2.before %>%
  add_column(sampleID = sampleID, .before = "Sample")
seq.r2.before <- seq.r2.before[,-2]

#seq.r2.after
seq.r2.after %>% as_tibble()
sampleID = gsub("(IGF[0-9]+)_.*", "\\1", seq.r2.after$Sample)
seq.r2.after <- seq.r2.after %>%
  add_column(sampleID = sampleID, .before = "Sample")

for (i in 1:length(seq.r2.after$Sample)){
  if (str_detect(seq.r2.after$Sample[i], "unpaired") == TRUE){
    seq.r2.after$Sample[i] <- "unpaired"
  }
  else
    seq.r2.after$Sample[i] <- "paired"
}


pivot_wider(seq.r2.after, names_from = Sample,
            values_from = FastQC_mqc.generalstats.fastqc.total_sequences)

#boxplot
ggplot(, aes(x=sample_type, y=segs, fill=sample_type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Number of segments in Primary Vs. Brain Mets") +
  xlab("")
