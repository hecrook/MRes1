library(tidyverse)
library(ggpubr)
library(rstatix)
#########RELATIVE PANEL - initial sequencing Vs. OCTOOPUS
# get_segments! #
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/Sequences_fig/")
list.files(getwd())
source("get_segments.R")
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/coverage_table/")
copyNumbersSegmented
x <- copyNumbersSegmented
segments.bm.1 <- get_segments(x)
segments.bm.1
write.table(segments.bm.1, file = "segments_BM_run1.txt")

#stratify segments for primary vs brain
igf.brain <- sid.IGF.anon$Sample_ID[which(sid.IGF.anon$sample_type == "brain")]
segments.bm.1$sample <- gsub("_sorted", "", segments.bm.1$sample)
segments.brain.r1 <- segments.bm.1 %>%
  filter(sample %in% igf.brain)
igf.primary <- sid.IGF.anon$Sample_ID[which(sid.IGF.anon$sample_type == "primary")]
segments.primary.r1 <- segments.bm.1 %>%
  filter(sample %in% igf.primary)
colnames(segments.primary.r1)[1] <- "sample.id"
colnames(segments.brain.r1)[1] <- "sample.id"

#get the error and the reads from copy numberscalled object
#expected error
View(copyNumbersCalled)
var <- expectedVariance(copyNumbersCalled) %>% sqrt
sample.id <-copyNumbersCalled@phenoData@data[["name"]]
error <- data.frame(var = var, sample.id = sample.id)
error <- as_tibble(error)
error
error$sample.id <- gsub("_sorted", "", error$sample.id)

# var <- expectedVariance(copyNumbersCalled) %>% sqrt
# sample.id <-copyNumbersCalled@phenoData@data[["name"]]
# error.obs <- data.frame(var = var, sample.id = sample.id)
# error.obs <- as_tibble(error.obs)
# error.obs
# error.obs$sample.id <- gsub("_sorted", "", error.obs$sample.id)

#get observed error
copynumber <- assayDataElement(x, "copynumber")[condition, , drop=FALSE]
sdFUN="sdDiffTrim"
sdFUN <- match.fun(sdFUN)
noise <- apply(copynumber, MARGIN=2L, FUN=sdFUN, na.rm=TRUE)
#add to error tibble
error <- error %>%
  add_column(error.obs = noise)
#rname the var column to expected error
colnames(error)[3] <- "error.obs"
colnames(error)[1] <- "error.exp"

error.brain <- error %>%
  filter(sample.id %in% igf.brain)
error.primary <- error %>%
  filter(sample.id %in% igf.primary)

reads <- copyNumbersCalled@phenoData@data[["total.reads"]]
reads <- data.frame(reads = reads, sample.id = gsub("_sorted", "", sample.id))
reads <- as_tibble(reads)
reads
reads$sample.id <- gsub("_sorted", "", reads$sample.id)
reads.brain <- reads %>%
  filter(sample.id %in% igf.brain)
reads.primary <- reads %>%
  filter(sample.id %in% igf.primary)

#merge all reads, error, and segs together into one table by sample ID
final.table.brain.run1 <- inner_join(segments.brain.r1, error.brain, by = "sample.id") %>%
  inner_join(., reads.brain, by = "sample.id")
final.table.primary.run1 <- inner_join(segments.primary.r1, error.primary, by = "sample.id") %>%
  inner_join(., reads.primary, by = "sample.id")

#t.test
# final.table.brain.run1 %>% summary(error.exp, type = "mean_sd")
# final.table.primary.run1 %>% summary(error.exp, type = "mean_sd")
y <- final.table.brain.run1$error.exp
z <- final.table.primary.run1$error.exp
t.test(y, z)
# Welch Two Sample t-test
# 
# data:  y and z
# t = -0.66174, df = 23.586, p-value = 0.5146
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.04622647  0.02379624
# sample estimates:
#   mean of x mean of y 
# 0.3795369 0.3907521 

#save outputs
write.csv(final.table.brain.run1, "Brain_Table_run1")
write.csv(final.table.primary.run1, "Primary_Table_run1")

#merge into 1 summary table for initial sequencing
final.table.run1 <- rbind(x = final.table.primary.run1, y = final.table.brain.run1)
write.csv(final.table.run1, "final_summary_run1")

##############topupseq summary
# need to pull out copynumberscalled and copynumberssegmented QDNASeq objects
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/Sequences_fig/")
source("get_segments.R")
# setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/coverage_table/summary_topupseq/data/")
# load("copyNumbersSegmented.RData")
# copyNumbersSegmented
# x <- copyNumbersSegmented
# segments.bm.2 <- get_segments(x)
# segments.bm.2
# write.table(segments.bm.2, file = "segments_BM_run2.txt")

# #stratify segments for primary vs brain
# igf.brain <- sid.IGF.anon$Sample_ID[which(sid.IGF.anon$sample_type == "brain")]
# segments.bm.2$sample <- gsub("_sorted", "", segments.bm.2$sample)
# segments.brain.r2 <- segments.bm.2 %>%
#   filter(sample %in% igf.brain)
# igf.primary <- sid.IGF.anon$Sample_ID[which(sid.IGF.anon$sample_type == "primary")]
# segments.primary.r2 <- segments.bm.2 %>%
#   filter(sample %in% igf.primary)
# colnames(segments.primary.r2)[1] <- "sample.id"
# colnames(segments.brain.r2)[1] <- "sample.id"

#pull out errors
# load("copyNumbersCalled.RData")
# View(copyNumbersCalled)
# var <- expectedVariance(copyNumbersCalled) %>% sqrt
# sample.id <-copyNumbersCalled@phenoData@data[["name"]]
# error <- data.frame(var = var, sample.id = sample.id)
# error <- as_tibble(error)
# error
# error$sample.id <- gsub("_sorted", "", error$sample.id)

# condition <- binsToUse(x)
# x <- copyNumbersCalled
# copynumber <- assayDataElement(x, "copynumber")[condition, , drop=FALSE]
# sdFUN="sdDiffTrim"
# sdFUN <- match.fun(sdFUN)
# noise <- apply(copynumber, MARGIN=2L, FUN=sdFUN, na.rm=TRUE)

############something wrong, incorrect segments being printed, going back to see if its the function 
# setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/coverage_table/summary_initialseq/")
# load("copyNumbersSegmented.RData")
# x <- copyNumbersSegmented
# segments.test <- get_segments(copyNumbersSegmented)
# that works abslutely fine, reloading the other data
# getting copynumberssegmented from top up sequencing, 30kb bins.
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/coverage_table/summary_topupseq/data/")
load("copyNumbersSegmented.RData")
x <- copyNumbersSegmented
segments.test.2 <- get_segments(copyNumbersSegmented)
#works absolutely fine, must have done something wrong previously, commented it out
# now go get the copynumberscalled from the same folder and save segments file
segments.test.2$sample <- gsub("_sorted", "", segments.test.2$sample)
colnames(segments.test.2)[1] <- "sample.id"
segments.brain.r2 <- segments.test.2 %>%
  filter(sample.id %in% igf.brain)
segments.primary.r2 <- segments.test.2 %>%
  filter(sample.id %in% igf.primary)

write.csv(segments.test.2, file = "segments.BM.run2.txt")
write.csv(segments.brain.r2, file = "segments.brain.run2")
write.csv(segments.primary.r2, file = "segments.primary.run2")

getwd()
load("copyNumbersCalled.RData")
x <- copyNumbersCalled


# observed error
copynumber <- assayDataElement(x, "copynumber")[condition, , drop=FALSE]
sdFUN="sdDiffTrim"
sdFUN <- match.fun(sdFUN)
noise <- apply(copynumber, MARGIN=2L, FUN=sdFUN, na.rm=TRUE)
noise

#expected error
var <- expectedVariance(copyNumbersCalled) %>% sqrt
sample.id <-copyNumbersCalled@phenoData@data[["name"]]
error.exp <- data.frame(var = var, sample.id = sample.id)
error.exp <- as_tibble(error.exp)
error.exp
error.exp$sample.id <- gsub("_sorted", "", error.exp$sample.id)
colnames(error.exp)[1] <- "error.exp"

#add observed error to error.exp
error.exp <- error.exp %>%
  add_column(error.obs = noise)
error <- as_tibble(error.exp)

#stratify into brain vs primary
error.brain <- error %>%
  filter(sample.id %in% igf.brain)
error.primary <- error %>%
  filter(sample.id %in% igf.primary)

reads <- copyNumbersCalled@phenoData@data[["total.reads"]]
reads <- data.frame(reads = reads, sample.id = gsub("_sorted", "", sample.id))
reads <- as_tibble(reads)
reads
reads$sample.id <- gsub("_sorted", "", reads$sample.id)
reads.brain <- reads %>%
  filter(sample.id %in% igf.brain)
reads.primary <- reads %>%
  filter(sample.id %in% igf.primary)

#merge all reads, error, and segs together into one table by sample ID
final.table.brain.run2 <- inner_join(segments.brain.r2, error.brain, by = "sample.id") %>%
  inner_join(., reads.brain, by = "sample.id")
final.table.primary.run2 <- inner_join(segments.primary.r2, error.primary, by = "sample.id") %>%
  inner_join(., reads.primary, by = "sample.id")

#save outputs
write.csv(final.table.brain.run2, "Brain_Table_run2")
write.csv(final.table.primary.run2, "Primary_Table_run2")

#merge into 1 summary table for initial sequencing
final.table.run2 <- rbind(x = final.table.primary.run2, y = final.table.brain.run2)
write.csv(final.table.run2, "final_summary_run2")

##################AGAIN for OCTOPUS
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/coverage_table/summary_OCTOPUS/")
load("copyNumbersSegmented.RData")
x <- copyNumbersSegmented
# segments.octopus <- get_segments(copyNumbersSegmented)
# only need archival samples so subset
segments.octopus$sample <- gsub("_sorted", "", segments.octopus$sample)
colnames(segments.octopus)[1] <- "sample.id"
# oct.meta <- read.csv("sWGS_metadata_HM.csv", header = T)
# oct.meta <- as_tibble(oct.meta)
# igf.arc <- oct.meta$IGF[which(oct.meta$SampleType == "Archival")]
#doesn't have as many rows as it should, should have 81 for all archival samples
oct.igf <- as_tibble(oct.igf)
oct.igf.2 <- as_tibble(oct.igf.2)
colnames(oct.igf)[3] <- "IMPCO"
oct.meta <- inner_join(oct.igf, oct.igf.2, by = "IMPCO")
oct.meta <- oct.meta[which(oct.meta$type.of.sample == "Arch"),]
igf.arc <- oct.meta$IGF[which(oct.meta$type.of.sample == "Arch")]

segments.arc <- segments.octopus %>%
  filter(sample.id %in% igf.arc)

write.csv(segments.octopus, file = "segments.oct.txt")
write.csv(segments.arc, file = "segments.arc")
write.csv(oct.meta, "OCTOPUS_meta_data")

#####errors
# need to load copynumberscalled object
getwd()
load("copyNumbersCalled.RData")
x <- copyNumbersCalled

# observed error
copynumber <- assayDataElement(x, "copynumber")[condition, , drop=FALSE]
sdFUN="sdDiffTrim"
sdFUN <- match.fun(sdFUN)
noise <- apply(copynumber, MARGIN=2L, FUN=sdFUN, na.rm=TRUE)
noise

#expected error
var <- expectedVariance(copyNumbersCalled) %>% sqrt
sample.id <-copyNumbersCalled@phenoData@data[["name"]]
error.exp <- data.frame(var = var, sample.id = sample.id)
error.exp <- as_tibble(error.exp)
error.exp
error.exp$sample.id <- gsub("_sorted", "", error.exp$sample.id)
colnames(error.exp)[1] <- "error.exp"

#add observed error to error.exp
error.exp <- error.exp %>%
  add_column(error.obs = noise)
error <- as_tibble(error.exp)

# #stratify into brain vs primary
# error.brain <- error %>%
#   filter(sample.id %in% igf.brain)
# error.primary <- error %>%
#   filter(sample.id %in% igf.primary)

reads <- copyNumbersCalled@phenoData@data[["total.reads"]]
reads <- data.frame(reads = reads, sample.id = gsub("_sorted", "", sample.id))
reads <- as_tibble(reads)
reads
reads$sample.id <- gsub("_sorted", "", reads$sample.id)
# reads.brain <- reads %>%
#   filter(sample.id %in% igf.brain)
# reads.primary <- reads %>%
#   filter(sample.id %in% igf.primary)

#merge all reads, error, and segs together into one table by sample ID
final.table.oct <- inner_join(segments.arc, error, by = "sample.id") %>%
  inner_join(., reads, by = "sample.id")
# final.table.primary.run2 <- inner_join(segments.primary.r2, error.primary, by = "sample.id") %>%
#   inner_join(., reads.primary, by = "sample.id")

#save outputs
write.csv(final.table.oct, "final_table_OCTOPUS")

# #merge into 1 summary table for initial sequencing
# final.table.run2 <- rbind(x = final.table.primary.run2, y = final.table.brain.run2)
# write.csv(final.table.run2, "final_summary_run2")

##merge all summary tables together
length(final.table.run1$sample.id) + length(final.table.run2$sample.id) + length(final.table.oct$sample.id)
# 132, this is how long the final table should be
# need to add_column to say which dataset they are coming from
premerge.table.run1 <- final.table.run1 %>%
  add_column(data = "Run1")
premerge.table.run2 <- final.table.run2 %>%
  add_column(data = "Run2")
premerge.table.oct <- final.table.oct %>%
  add_column(data = "OCT.")
#merge with rbind()
final.table <- rbind(premerge.table.run1, premerge.table.run2, premerge.table.oct)
write.csv(final.table, file = "final_table.csv")

#calculate coverage
final.table$coverage <- (final.table$reads *150)/3000000000
write.csv(final.table, file = "final_table_coverage.csv")

#create mean of octopus data for coverage table
final.table.oct$coverage <- (final.table.oct$reads * 150)/3000000000
means.final.table.oct <- data.frame(n.segments = NA, error.exp = NA, error.obs = NA, reads = NA, coverage = NA)
means.final.table.oct <- means.final.table.oct %>%
  means.final.table.oct$n_segments <- mean(final.table.oct$n_segments)
  means.final.table.oct$error.exp <- mean(final.table.oct$error.exp) +
  means.final.table.oct$error.obs <- mean(final.table.oct$error.obs) +
  means.final.table.oct$reads <- mean(final.table.oct$reads) +
  means.final.table.oct$coverage <- mean(final.table.oct$coverage)

visual <- final.table.run1
visual$coverage %>% mean(visual$reads)
visual$coverage <- mean(visual$reads)
rm(visual$coverage)
visual
visual <- visual[,-6]
visual
visual$coverage <- (visual$reads * 150)/3000000000
visual
library(tidyverse)
visual <- visual %>%
visual$reads / 1000000
visual$reads / 1000000
visual$reads <- visual$reads / 1000000
View(visual)
visual <- visual[,c(5,6,2,3,4)]
write.csv(visual, file = "visual.csv")
getwd()
round(visual$reads, 2)
visual$reads <- round(visual$reads, 2)
visual$coverage <- round(visual$coverage, 2)
visual$error.exp <- round(visual$error.exp, 2)
visual$error.obs <- round(visual$error.obs, 2)
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/coverage_table/")
write.csv(visual, file = "visual.csv")
# plot with mmtable2
# remotes::install_github("ianmoran11/mmtable2")
# library(mmtable2)
# library(gt)
# library(tidyverse)
# 
# table <- final.table %>%
#   mmtable(table_data = reads, table_name = "reads")

#boxlpot of coverage 
# ggplot(Brain.gg, aes(x=sample_type, y=segs, fill=sample_type)) +
#   geom_boxplot() +
#   scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
#   geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle("Number of segments in Primary Vs. Brain Mets") +
#   xlab("")

library(tidyverse)
library(hrbrthemes)
library(viridis)
boxplot <- final.table[,c(6:7)]
ggplot(boxplot, aes(x=data, y=coverage, fill=data)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  ggtitle("coverage") +
  xlab("")

write.csv(final.table, file = "final_table_wcoverage.csv")
mean(boxplot$coverage[which(boxplot$data == "Run2")])
#looks bad, try bar graph
ggplot(boxplot) +
  geom_bar( aes(x=data, y=coverage), stat = "identity", fill = "skyblue", alpha = 0.7) +
  scale_x_discrete(limits = c("OCT.", "Run1", "Run2"))

#try just the two different runs
run1vrun2 <- rbind(premerge.table.run1, premerge.table.run2)
run1vrun2$coverage <- (run1vrun2$reads *150)/3000000000
boxplot <- run1vrun2[,c(6:7)]
ggplot(boxplot, aes(x=data, y=coverage, fill=data)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  ggtitle("coverage") +
  xlab("")

#try sequences

boxplot <- final.table[,c(4,6)]
ggplot(boxplot, aes(x=data, y=error.obs, fill=data)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  ggtitle("Observed error across experiments") +
  xlab("")
