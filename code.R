library(QDNAseq)
library(Biobase)
library(ACE)
library(dplyr)
library(GenomicRanges)

setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/ACE/")
load("copyNumbersCalled.RData")

#Processing CN metadata
file1 <- read.table("mcneish_2-11-2020_sWGS_edited_HM.txt", header = T, sep = "\t", stringsAsFactors = F)
file2 <- read.table("sWGS_PatientID_Type_IMPCO_HM.txt", header = T, sep = "\t", stringsAsFactors = F)
#changes "anonymized dolcumn" header to IMPCO so that we can merge file 1 and file 2 based on the IMPCO column matches
colnames(file2)[3] <- "IMPCO"
meta <- merge(file1,file2,by="IMPCO", all = T)
rm(file1,file2)
#Changing column names of meta data so it can be merged with tp53 data
colnames(meta)[3:4] <- c("PatientID", "SampleType")
meta$SampleType[which(meta$SampleType=="Arch")] <- "Archival"
meta$SampleType[which(meta$SampleType=="SE")] <- "StudyEntry"
meta$SampleType[which(meta$SampleType=="wb")] <- "WholeBlood"

#Loading TP53 VAF
tp53 <- read.table("Final_TP53_VAF_HM_v2.txt", header = T, sep = "\t", stringsAsFactors = F)
tp53 <- merge(meta, tp53, by=c("PatientID", "SampleType"), all = T)
#Remove instances with missing IMPCO, IGF, or PatientID
tp53 <- tp53[-which(is.na(tp53$IMPCO)),]
tp53 <- tp53[-which(is.na(tp53$PatientID)),]
tp53 <- tp53[-which(is.na(tp53$IGF)),]
#rm any samples where there is no available tp53 VAF from either Mutect2 AND Strelka
tp53.avail <- tp53[-intersect(which(is.na(tp53$Mutect2)),which(is.na(tp53$Strelka))),]
#create a column in tp53.avail with the difference between the two VAF scores
tp53.avail$vaf.diff <- abs(tp53.avail$Mutect2 - tp53.avail$Strelka)
#Remove the instances where there are IGF samples not in copyNumbersCalled
tp53.avail.OHA <- tp53.avail[tp53.avail$IGF %in% c("IGF116720", "IGF116728", "IGF116737", "IGF116740", "IGF116745", "IGF116760", "IGF116768", "IGF116776", "IGF116779", "IGF116798", "IGF116816", "IGF116833", "IGF116836", "IGF116840", "IGF116842", "IGF116844", "IGF116849", "IGF116854", "IGF116856", "IGF116865"),]

#Fixing the CN object. sampleNames() access the feature names and sample names stored in an object.
sampleNames(copyNumbersCalled) <- sapply(strsplit(as.character(sampleNames(copyNumbersCalled)), "_"), "[[", 1)
pData(copyNumbersCalled)$name <- sapply(strsplit(as.character(pData(copyNumbersCalled)$name), "_"), "[[", 1)
copyNumbersCalled.avail <- copyNumbersCalled[,tp53.avail.OHA$IGF]
#Calculation to check the data are the same.
all(sampleNames(copyNumbersCalled.avail)==tp53.avail.OHA$IGF)
all(pData(copyNumbersCalled.avail)$name==tp53.avail.OHA$IGF)
#merge columns from both objects
pData(copyNumbersCalled.avail) <- cbind(pData(copyNumbersCalled.avail),tp53.avail.OHA)


# ACE using squaremodel fitting (code adopted from Poell et al i.e. original ACE paper)
# segrel <- readRDS("OVCAR_rel.rds")
segrel <- copyNumbersCalled.avail
# # Manual calls!!
#create a data frame which will contain info on model picked for all samples
#manual.fits <- data.frame("ploidy"=NA,"cellularity"=NA,"error"=NA,"minimum"=NA)
#manual.fits <- manual.fits[-1,]
# # To run for each sample manually, run the line below by changing index
s=20
cat("Working with ", sampleNames(segrel)[s]," | ",s,"/",length(sampleNames(segrel)),"\n")
minimadf <- squaremodel(segrel, s, ptop = 4.3, pbottom = 1.8, prows = 250, penalty = 0.5, penploidy = 0.5)$minimadf
View(minimadf)

#view cellularity of sample from tp53 VAF 
sampleNames(segrel)[s]
tp53.avail.OHA[tp53.avail.OHA$IGF == 'IGF116856',]

#select model and put into manual.fits data frame
selected.row <- minimadf['22105',]
manual.fits <- rbind(manual.fits,selected.row)
rownames(manual.fits)[nrow(manual.fits)] <- sampleNames(segrel)[s]

to.rm <- c("IGF116775","IGF116726","IGF116784","IGF116738","IGF116809","IGF116841","IGF116827","IGF116858","IGF116746","IGF116835","IGF116860","IGF116805","IGF116807","IGF116800")
copyNumbersCalled.final <- segrel[,setdiff(sampleNames(segrel),to.rm)]

#save
save(manual.fits,file="manual.fits.RData")


###########################################################################

#Writing/saving outputs
dir.create("Final_CN_OHA")
write.table(manual.fits,file="Final_CN_OHA/manual.fits.txt", row.names = T, col.names = NA, sep ="\t", quote = F)
save(manual.fits, file = "Final_CN_OHA/manual.fits.RData")

write.table(tp53.avail, file="Final_CN_OHA/tp53.avail.txt", row.names = F, col.names = T, sep = "\t", quote = F)
saveRDS(copyNumbersCalled.avail, file = "Final_CN_OHA/copynumbersCalled.avail.rds")

pData(copyNumbersCalled.avail)$name <- gsub("_sorted","",pData(copyNumbersCalled.avail)$name)
saveRDS(copyNumbersCalled.avail, file = "Final_CN_OHA/copyNumbersCalled.avail.rds")
models.file <- data.frame("Sample"=rownames(manual.fits),"Cellularity"=manual.fits$cellularity,"Ploidy"=manual.fits$ploidy)
write.table(models.file, "Final_CN_OHA/models.file.tsv", row.names = FALSE, sep = "\t", quote = FALSE)

##########################################################################
# Fitting ACE based on above variant
# had to remove outputdir argument eventually >:(
# Also was orignially copynumbersCalled.final.rds which had the samples that were demmed unsueable by Hasans standards (something to do with tp52 VAF), but that was causing issues with NAs - probably something to do with the data frames not matching up perfectly
postanalysisloop("Final_CN_OHA/copynumbersCalled.avail.rds", modelsfile = "Final_CN_OHA/models.file.tsv", trncname = FALSE, ouputdir = "Final_CN_OHA/ACE_output")

#checks something? idk what but its true so we move
all(names(copyNumbersCalled.avail)==models.file$Sample)

###########################################################################
### Input for CN signatures
setwd("~/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/ACE/segmentfiles/")
#CN signature input files are the segmentfiles 
files <- list.files(getwd(),pattern = "_segments.tsv")
samples <- gsub("_segments.tsv","",files)
#tab.lis becomes a list of dataframes for each CN
tab.lis <- list()

########FOR LOOP########
for(i in 1:length(files)) {
  tab <- read.table(files[i], header = T, sep = "\t", stringsAsFactors = F)
  tab.lis[[i]] <-tab[,c(1:3,5)]
  colnames(tab.lis[[i]]) <- c("chromosome", "start", "end", "segVal")
}

names(tab.lis) <- samples
save(tab.lis, file = "tab.lis.ace.RData")

# CN signatures


##############################PAUSE, BACK TO SCRIPT ONE TO CREATE FINAL_TABLE_EVERYTHING.txt
#merge models.fit and tp53 avail
models.file <- read.table("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/ACE/Final_CN_OHA/models.file.tsv", header = T, sep = "\t", stringsAsFactors = F)
tp53.avail <- read.table("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/ACE/Final_CN_OHA/tp53.avail.txt", header = T, sep = "\t", stringsAsFactors = F)
colnames(models.file)[2:3] <- paste("ACE_", colnames(models.file)[2:3],sep = "")
colnames(models.file)[1] <- "IGF"
tp53.avail <- merge(tp53.avail,models.file,by="IGF")

final.table.everything.OHA <- tp53.avail
write.table(final.table.everything.OHA, file = "final_table_everything_OHA.txt", row.names = F, col.names = T, sep = "\t", quote = F)
save(final.table.everything.OHA, file = "final_table_everything_OHA.RData")
#################################BACK TO SCRIPT 2

load("final_table_everything_OHA.RData")
for(i in 1:ncol(octopus.cat)) {
  cat("Working with ",colnames(octopus.cat)[i], " | ",i,"/",ncol(octopus.cat),"\n")
  pl <- final.table.everything.OHA$ACE_Ploidy[which(final.table.everything.OHA$IGF==colnames(octopus.cat)[i])]
  if(pl<=2.7) {
    octopus.cat[,i]  <- "Normal"
    octopus.cat[which(temp[,i]<1),i]  <- "Loss"
    octopus.cat[which(temp[,i]>=5),i]  <- "Gain"
    
  } else if(pl>2.7) {
    magicnumber <- pl-2.7
    octopus.cat[,i]  <- "Normal"
    octopus.cat[which(temp[,i]<magicnumber),i]  <- "Loss"
    octopus.cat[which(temp[,i]>=9),i]  <- "Gain"
  }
}

# Calculating frequencies (PERCENTAGES)
freq.loss <- freq.normal <- freq.gain <- NULL
for(i in 1:nrow(octopus.cat)) {
  freq.loss <- c(freq.loss,round((length(which(octopus.cat[i,]=="Loss"))/ncol(octopus.cat))*100))
  freq.normal <- c(freq.normal,round((length(which(octopus.cat[i,]=="Normal"))/ncol(octopus.cat))*100))
  freq.gain <- c(freq.gain,round((length(which(octopus.cat[i,]=="Gain"))/ncol(octopus.cat))*100))
}
octopus.freq <- OHA.consensus[,1:4]
octopus.freq$Loss <- freq.loss
octopus.freq$Normal <- freq.normal
octopus.freq$Gain <- freq.gain

# Calculating frequencies (SAMPLE NUMBER)
freq.loss <- freq.normal <- freq.gain <- NULL
for(i in 1:nrow(octopus.cat)) {
  freq.loss <- c(freq.loss,length(which(octopus.cat[i,]=="Loss")))
  freq.normal <- c(freq.normal,length(which(octopus.cat[i,]=="Normal")))
  freq.gain <- c(freq.gain,length(which(octopus.cat[i,]=="Gain")))
}
octopus.freq.num <- OHA.consensus[,1:4]
octopus.freq.num$Loss <- freq.loss
octopus.freq.num$Normal <- freq.normal
octopus.freq.num$Gain <- freq.gain

####################################
## Genome plots - Frequency plots ##
####################################

# colo <- c("firebrick3","darkgoldenrod1","chartreuse4","deepskyblue","darkorange1","darkorchid","lemonchiffon4","darkgrey","magenta")
colo <- c("firebrick3","cornflowerblue","darkblue")
chr <- c(1:22,"X")
CNtypes <- c("Loss","Gain")

library(grDevices)

#Cairo PDF
GroupName <- "OHA"
freq.tab <- octopus.freq[,-which(colnames(octopus.freq)=="Normal")]
cairo_pdf(paste("FrequencyPlots_All_CNtypes_",GroupName,"_v2_COSMIC.pdf",sep=""), width = 14, height = 10, onefile=TRUE)
par(mfrow=c(2,1))
for(j in 1:length(CNtypes)) {
  plot(x=xlimit,y=c(0,100),type="n", xaxt="n", yaxt="n",xlab = "", ylab=paste(CNtypes[j]," Percentage",sep=""),main = paste("OCTOPUS ",CNtypes[j],sep=""))
  axis(side=2, at=seq(0,100,by=25), cex.axis=1,las=1,labels=c(0,25,50,75,100))
  axis(side=1, at=meanind, cex.axis=0.7,las=1,labels=chr,lwd.ticks=1)
  for (i in 1:nrow(freq.tab)) {
    if(freq.tab$chromosome[i]=="1") {
      rect(freq.tab$start[i],0,freq.tab$end[i],freq.tab[i,4+j],col=colo[j],border=NA)
    }
    else{
      start <- freq.tab$start[i]+ref$LengthAdded[which(ref$Chromosome==freq.tab$chromosome[i])-1]
      end <- freq.tab$end[i]+ref$LengthAdded[which(ref$Chromosome==freq.tab$chromosome[i])-1]
      rect(start,0,end,freq.tab[i,4+j],col=colo[j],border=NA)
    }
  }
  # abline(v=allind,lty="dotted",lwd=0.3,a=1,b=2)
  abline(h=0,lty="solid",lwd=1)
  abline(h=c(-100,-75,-50,-25,0,25,50,75,100),lty="dotted",lwd=0.3)
  abline(v=ref$End,lty="dotted",lwd=0.3)
}
dev.off()
