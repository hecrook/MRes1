setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/MultiQC results/")
multiqc <- as_tibble(read.delim("multiqc_general_stats.txt"))
multiqc
colnames(multiqc)
mean.seq <- mean(multiqc$FastQC_mqc.generalstats.fastqc.total_sequences)
mean.seq
# multiqc.fastqc <- as_tibble(read.delim("multiqc_fastqc.txt"))
# mean.adap <- mean(multiqc.fastqc$adapter_content)

setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/MultiQC results/before/")
multi.after <- as_tibble(read.delim("multiqc_general_stats.txt"))
mean.seq <- mean(multi.after$FastQC_mqc.generalstats.fastqc.total_sequences)
(mean.seq * 150)/3000000000
(1317334 * 150)/3000000000
mean.after.r1