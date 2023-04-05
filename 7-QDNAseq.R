# Hannah Crook
# 09/03/2023

# Running QDNAseq on sWGS OHA

# load required packages
library(QDNAseq)
library(Biobase)
library(DNAcopy)
library(CGHcall)
library(future)

# Set working directory
setwd("/rds/general/user/hec22/projects/mcneish_team_data/ephemeral/IGFQ001515_mcneish_14-12-2022_lcWGS/data/reads2/")

# To process data in parallel using multiple processes on the current machine, use the following. After that, all functions that support parallel processing will automatically use it.
future::plan("multisession", workers=4)

# Then we need to obtain bin annotations. These are available pre-calculated for genome build
# mm10 and bin sizes 1, 5, 10, 15, 30, 50, 100, 500, and 1000 kbp. 
bins <- getBinAnnotations(binSize=30,genome="hg19")

### Main loop (per sample analysis)
samples <- list.dirs(getwd(),recursive=F,full.names=F)
bamfiles <- list.files(getwd(),pattern = "_sorted.bam",recursive=T)
bamfiles <- bamfiles[-grep(".bai",bamfiles)]
cat("Reading counts per bin ..\n")
readCounts <- binReadCounts(bins = bins,bamfiles = bamfiles, cache = T,pairedEnds = F,isDuplicate=F)
save(readCounts,file="readCounts.RData")

# Plot a raw copy number profile (read counts across the genome), and highlight bins that will be removed with default filtering
cat("Plotting raw copy number profiles i.e. read counts across the genome..\n")
cairo_pdf("RawCopyNumberHighlights.pdf",width = 10, height = 7,onefile = T)
plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
dev.off()

# Apply filters and plot median read counts as a function of GC content and mappability
cat("Applying filters and plotting median read counts as a function of GC content and mappability..\n")
readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
cairo_pdf("isobarPlot.pdf",width = 10, height = 7,onefile = T)
isobarPlot(readCountsFiltered)
dev.off()
save(readCountsFiltered,file="readCountsFiltered.RData")

# Estimate the correction for GC content and mappability, and make a plot for the relationship between the observed standard deviation in the data and its read depth. The theoretical expectation is a linear relationship, which is shown in the plot with a black line. Samples with low-quality DNA will be noisier than expected and appear further above the line than good-quality samples
cat("Plotting observed SD vs read depth..\n")
readCountsFiltered <- estimateCorrection(readCountsFiltered)
save(readCountsFiltered,file="readCountsFiltered_correctionEstimated.RData")
cairo_pdf("noisePlot.pdf",width = 10, height = 7,onefile = T)
noisePlot(readCountsFiltered)
dev.off()

# Next, we apply the correction for GC content and mappability. This will return a QDNAseqCopyNumbers object, which we then normalize, smooth outliers, and plot the copy number profile
cat("GC corrections..\n")
copyNumbers <- correctBins(readCountsFiltered)
save(copyNumbers,file="copyNumbers.RData")
cat("Normalisation..\n")
copyNumbersNormalized <- normalizeBins(copyNumbers)
save(copyNumbersNormalized,file="copyNumbersNormalized.RData")
cat("Smoothing outliers..\n")
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
save(copyNumbersSmooth,file="copyNumbersSmooth.RData")
cat("Plotting smooothened CN profiles..\n")
cairo_pdf("copyNumbersSmooth.pdf",width = 10, height = 7,onefile = T)
plot(copyNumbersSmooth)
dev.off()

# Data is now ready to be analyzed with a downstream package of choice. For analysis with an external program or for visualizations in IGV, the data can be exported to a file
cat("Exporting binned data in three different formats..\n")
exportBins(copyNumbersSmooth, file="copyNumbersSmooth.txt")
exportBins(copyNumbersSmooth, file="copyNumbersSmooth.igv", format="igv")
#exportBins(copyNumbersSmooth, file="copyNumbersSmooth.bed", format="bed")

# Segmentation with the CBS algorithm from DNAcopy, and calling copy number aberrations with CGHcall or cutoffs have been implemented for convenience. By default, segmentation uses a log2-transformation, but a sqrt(x + 3/8) can also be used as it stabilizes the variance of a Poisson distribution (Anscombe transform)
cat("CBS segmentation..\n")
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
save(copyNumbersSegmented,file="copyNumbersSegmented.RData")
cat("Normalising CBS segmented data..\n")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
save(copyNumbersSegmented,file="copyNumbersSegmented_normalised.RData")
cat("Plotting normalised CBS segmented data..\n")
cairo_pdf("copyNumbersSegmented.pdf",width = 10, height = 7,onefile = T)
plot(copyNumbersSegmented)
dev.off()
cat("Calling CN per bins..\n")
copyNumbersCalled <- callBins(copyNumbersSegmented)
save(copyNumbersCalled,file="copyNumbersCalled.RData")
cairo_pdf("copyNumbersCalled.pdf",width = 10, height = 7,onefile = T)
plot(copyNumbersCalled)
dev.off()

# HM
# Saving segment level CN status probabilities
copyNumbersCalled_calls <- assayData(copyNumbersCalled)$calls
save(copyNumbersCalled_calls,file="copyNumbersCalled_calls.RData")
copyNumbersCalled_copynumber <-assayData(copyNumbersCalled)$copynumber
save(copyNumbersCalled_copynumber,file="copyNumbersCalled_copynumber.RData")
copyNumbersCalled_probamp <-assayData(copyNumbersCalled)$probamp
save(copyNumbersCalled_probamp,file="copyNumbersCalled_probamp.RData")
copyNumbersCalled_probdloss <-assayData(copyNumbersCalled)$probdloss
save(copyNumbersCalled_probdloss,file="copyNumbersCalled_probdloss.RData")
copyNumbersCalled_probgain <-assayData(copyNumbersCalled)$probgain
save(copyNumbersCalled_probgain,file="copyNumbersCalled_probgain.RData")
copyNumbersCalled_probloss <-assayData(copyNumbersCalled)$probloss
save(copyNumbersCalled_probloss,file="copyNumbersCalled_probloss.RData")
copyNumbersCalled_probnorm <-assayData(copyNumbersCalled)$probnorm
save(copyNumbersCalled_probnorm,file="copyNumbersCalled_probnorm.RData")
copyNumbersCalled_segmented <-assayData(copyNumbersCalled)$segmented
save(copyNumbersCalled_segmented,file="copyNumbersCalled_segmented.RData")
