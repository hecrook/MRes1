#Creating finaltableeverything
library(stringr)
final <- sid.IGF.anon
colnames(final) <- c("ANON", "Gaia_cellularity", "Sample", "Sampletype")
final <- inner_join(final, models.file, by = "Sample")
# Create a patient ID column
final$PatientID <- NA
for (i in 1:length(final$ANON)){
  if (sum(str_detect(final$ANON[i], 'A')) >0){
    final[i,7] <- gsub("A", "", as.character(final[i,1]))
  } else {
    final[i,7] <- gsub("B", "", as.character(final[i,1]))
  }
}
# Copy Number GAINS AND LOSSES
library(GenomicRanges,quietly=T)
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/BrainMets_CNS2/data/100kb/segmentfiles/")
files <- list.files(getwd(),pattern = "_segments.tsv")
samples <- gsub("_segments.tsv","",files)
tab.lis <- list()
for(i in 1:length(files)) {
  tab <- read.table(files[i],header = T,sep = "\t",stringsAsFactors = F)
  tab.lis[[i]] <- tab[,c(1:3,5)]
  colnames(tab.lis[[i]]) <- c("chromosome", "start", "end", "segVal")
}

# lapply(tab.lis,head)
names(tab.lis) <- samples 
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/BrainMets_CNS2/data/100kb/")
save(tab.lis,file="tab.lis.ace.RData")

cat("Getting all breakpoints..\n")
Breakpoints <- vector("list", 22)
# chr <- c(1:22,"X")
chr <- c(1:22)
for (i in 1:length(tab.lis)) {
  cat("Working with ",names(tab.lis)[i], " | ",i,"/",length(tab.lis),"\n")
  file <- tab.lis[[i]]
  file$start <- as.numeric(file$start)
  file$end <- as.numeric(file$end)
  for(j in 1:22) {
    vec <- 0
    ind <- which(file$chromosome==chr[j])
    vec<- c(vec,file$start[ind],file$end[ind])
    Breakpoints[[j]] <- c(Breakpoints[[j]],vec)
  }
}

# Sorting all breakpoints per chromosome after removing redundant breakpoints
cat("Sorting breakpoints per chromosome and removing duplicates..\n")
for(j in 1:length(Breakpoints)) {
  Breakpoints[[j]] <- sort(unique(Breakpoints[[j]]))
}

# Creating a data frame and turning all breakpoints into start and end of segments i.e. creating consensus segments
cat("Creating consensus segments from breakpoints..\n")
vec <- 0
vec <- vec[-1]
for(i in 1:length(Breakpoints)) {
  vec <- c(vec,rep(chr[i],length(Breakpoints[[i]])))
}
newcolnames <- names(tab.lis)
mat <- matrix(nrow=length(vec),ncol=length(newcolnames)+3)
df <- as.data.frame(mat)
colnames(df) <- c("chromosome","start","end",newcolnames)
df$chromosome <- vec
tag <- 0
for(i in 1:length(Breakpoints)) {
  # cat(i,"\n")
  bk <- Breakpoints[[i]]
  for(j in 1:(length(Breakpoints[[i]])-1))
  {
    # cat(j,"\n")
    df$start[tag+j] <- bk[j]
    df$end[tag+j] <- bk[j+1]
  }
  tag <- tag+length(Breakpoints[[i]])
}
df <- df[complete.cases(df$start),]
# Mapping CN from all samples to consensus segments (Code Corrected with Jelmar 21-07-14)
cat("Mapping CN from all TAPS CN files to consensus segments..\n")
# consensusG <- GRanges(seqnames = gsub("chr", "", df$chromosome, fixed=TRUE),ranges = IRanges(df$start, df$end))
consensusG <- GRanges(seqnames = df$chromosome,ranges = IRanges(df$start, df$end))
count <- 0
df2 <- df[,1:3]
for(i in 1:length(tab.lis)) {
  cat("Working with ",names(tab.lis)[i], " | ",i,"/",length(tab.lis),"\n")
  file <- tab.lis[[i]]
  file$start <- as.numeric(file$start)
  file$end <- as.numeric(file$end)
  colnames(file)[ncol(file)] <- "aCN"
  file$aCN <- as.character(file$aCN)
  fileG <- GRanges(seqnames = file$chromosome,ranges = IRanges(file$start, file$end), Cn=file$aCN)
  temp <- findOverlaps(consensusG, fileG, type="within")
  
  df2$Cn[queryHits(temp)] <- file$aCN[subjectHits(temp)]
  colnames(df2)[ncol(df2)] <- names(tab.lis)[i]
}

# Adding length and column order/names
var <- 0
for(i in 1:nrow(df2)) {
  var[i] <- (df2$end[i]-df2$start[i])/10^6
}
df2$LengthMB <- var
df2 <- df2[,c(1:3,ncol(df2),4:(ncol(df2)-1))]

# Removing segments with just 1 base (technical artifact of how I do consensus segments!)
df2 <- df2[-which(df2$LengthMB==0.000001),]

# First segment starting from 0 for each chromosome is usually NA, removing that
df2[which(df2$start==0),] # manually checked that all are NAs in that row
df2 <- df2[-which(df2$start==0),]

brainmets.consensus <- df2

#### CALCULATING FREQUENCIES PER GROUP
brainmets.cat <- as.data.frame(matrix(nrow = nrow(brainmets.consensus),ncol = ncol(brainmets.consensus)-4))
colnames(brainmets.cat) <- colnames(brainmets.consensus)[5:ncol(brainmets.consensus)]
# Rounding off to 2 digits
temp <- brainmets.consensus[,5:ncol(brainmets.consensus)]
for(i in 1:ncol(temp)) {
  temp[,i] <- round(as.numeric(temp[,i]),2)
}
# Upscaling and downscaling values to closest whole value if near than 10% 
for(i in 1:ncol(temp)) {
  temp[which(temp[,i]>=1.9 & temp[,i]<2),i] <- 2
  temp[which(temp[,i]>=2.9 & temp[,i]<3),i] <- 3
  temp[which(temp[,i]>=3.9 & temp[,i]<4),i] <- 4
}

# New definitions
# Categorising CN as per COSMIC guidelines
# https://cancer.sanger.ac.uk/cosmic/help/cnv/overview
# Ploidy <= 2.7
# Gain:
# average genome ploidy <= 2.7 AND total copy number >= 5
# Loss:
#   average genome ploidy <= 2.7 AND total copy number = 0
# Ploidy > 2.7
# Gain
# total copy number >= 9
# Loss
# total copy number < ( average genome ploidy - 2.7 )

#load("final_table_everything_HM.RData")
for(i in 1:ncol(brainmets.cat)) {
  cat("Working with ",colnames(brainmets.cat)[i], " | ",i,"/",ncol(brainmets.cat),"\n")
  pl <- final$ACE_Ploidy[which(final$IGF==colnames(brainmets.cat)[i])]
  if(pl<=2.7) {
    brainmets.cat[,i]  <- "Normal"
    brainmets.cat[which(temp[,i]<1),i]  <- "Loss"
    brainmets.cat[which(temp[,i]>=5),i]  <- "Gain"
    
  } else if(pl>2.7) {
    magicnumber <- pl-2.7
    brainmets.cat[,i]  <- "Normal"
    brainmets.cat[which(temp[,i]<magicnumber),i]  <- "Loss"
    brainmets.cat[which(temp[,i]>=9),i]  <- "Gain"
  }
}
# Calculating frequencies for whole cohort (PERCENTAGES)
freq.loss <- freq.normal <- freq.gain <- NULL
for(i in 1:nrow(brainmets.cat)) {
  freq.loss <- c(freq.loss,round((length(which(brainmets.cat[i,]=="Loss"))/ncol(brainmets.cat))*100))
  freq.normal <- c(freq.normal,round((length(which(brainmets.cat[i,]=="Normal"))/ncol(brainmets.cat))*100))
  freq.gain <- c(freq.gain,round((length(which(brainmets.cat[i,]=="Gain"))/ncol(brainmets.cat))*100))
}
brainmets.freq <- brainmets.consensus[,1:4]
brainmets.freq$Loss <- freq.loss
brainmets.freq$Normal <- freq.normal
brainmets.freq$Gain <- freq.gain

# Calculating frequencies (SAMPLE NUMBER)
freq.loss <- freq.normal <- freq.gain <- NULL
for(i in 1:nrow(brainmets.cat)) {
  freq.loss <- c(freq.loss,length(which(brainmets.cat[i,]=="Loss")))
  freq.normal <- c(freq.normal,length(which(brainmets.cat[i,]=="Normal")))
  freq.gain <- c(freq.gain,length(which(brainmets.cat[i,]=="Gain")))
}
brainmets.freq.num <- brainmets.consensus[,1:4]
brainmets.freq.num$Loss <- freq.loss
brainmets.freq.num$Normal <- freq.normal
brainmets.freq.num$Gain <- freq.gain

####################################
## Genome plots - Frequency plots ##
####################################
# old.par <- par()
load("hg19_stats.Rdata")
library(Cairo)
library(tibble)

# colo <- c("firebrick3","darkgoldenrod1","chartreuse4","deepskyblue","darkorange1","darkorchid","lemonchiffon4","darkgrey","magenta")
colo <- c("firebrick3","cornflowerblue","darkblue")
chr <- c(1:22,"X")
CNtypes <- c("Loss","Gain")

#stratify patients based on brainmets vs primary
primary.cat <- brainmets.cat %>%
  select(final %>% filter(SampleType == 'primary') %>% pull(IGF))
brain.cat <- brainmets.cat %>%
  select(final %>% filter(SampleType == 'brain') %>% pull(IGF))

# Cairo PDF
GroupName <- "BrainMets"
freq.tab <- brainmets.freq[,-which(colnames(brainmets.freq)=="Normal")]
cairo_pdf(paste("FrequencyPlots_All_CNtypes_",GroupName,"_v2_COSMIC.pdf",sep=""), width = 14, height = 10, onefile=TRUE)
par(mfrow=c(2,1))
for(j in 1:length(CNtypes)) {
  plot(x=xlimit,y=c(0,100),type="n", xaxt="n", yaxt="n",xlab = "", ylab=paste(CNtypes[j]," Percentage",sep=""),main = paste("BrainMets",CNtypes[j],sep=""))
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

# Merging together
GroupName <- "BrainMets"
freq.tab <- brainmets.freq[,-which(colnames(brainmets.freq)=="Normal")]
freq.tab$Loss <- freq.tab$Loss*-1
cairo_pdf(paste("FrequencyPlots_All_CNtypes_",GroupName,"_combined_v2_COSMIC.pdf",sep=""), width = 14, height = 5, onefile=TRUE)
par(mfrow=c(1,1))
plot(x=xlimit,y=c(-100,100),type="n", xaxt="n", yaxt="n",xlab = "", ylab="Percentage of CN change",main = "BrainMets, N = 24 | Frequency of copy number changes")
axis(side=2, at=seq(-100,100,by=25), cex.axis=1,las=1,labels=c(100,75,50,25,0,25,50,75,100))
axis(side=1, at=meanind, cex.axis=0.7,las=1,labels=chr,lwd.ticks=1)
for(j in 1:length(CNtypes)) {
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
}
# abline(v=allind,lty="dotted",lwd=0.3,a=1,b=2)
abline(h=0,lty="solid",lwd=1)
abline(h=c(-100,-75,-50,-25,0,25,50,75,100),lty="dotted",lwd=0.3)
abline(v=ref$End,lty="dotted",lwd=0.3)
legend("topright",legend = CNtypes,col = colo,lwd = 2,cex = 0.7)
dev.off()

################################################
# GROUP COMPARISONS | STATISTICAL SIGNIFICANCE #
################################################
library(stats)

load("final_table_everything_HC.RData")

primary <- final$IGF[which(final$SampleType=="primary")]
brain <- final$IGF[which(final$SampleType=="brain")]

# Calculating frequencies (SAMPLE NUMBER)
tempdf <- brainmets.cat[,which(colnames(brainmets.cat) %in% primary)]
freq.loss <- freq.normal <- freq.gain <- NULL
for(i in 1:nrow(tempdf)) {
  freq.loss <- c(freq.loss,length(which(tempdf[i,]=="Loss")))
  freq.normal <- c(freq.normal,length(which(tempdf[i,]=="Normal")))
  freq.gain <- c(freq.gain,length(which(tempdf[i,]=="Gain")))
}
brainmets.freq.num <- brainmets.consensus[,1:4]
brainmets.freq.num$Loss.primary <- freq.loss
brainmets.freq.num$Normal.primary <- freq.normal
brainmets.freq.num$Gain.primary <- freq.gain
tempdf <- brainmets.cat[,which(colnames(brainmets.cat) %in% brain)]
freq.loss <- freq.normal <- freq.gain <- NULL
for(i in 1:nrow(tempdf)) {
  freq.loss <- c(freq.loss,length(which(tempdf[i,]=="Loss")))
  freq.normal <- c(freq.normal,length(which(tempdf[i,]=="Normal")))
  freq.gain <- c(freq.gain,length(which(tempdf[i,]=="Gain")))
}
brainmets.freq.num$Loss.brain <- freq.loss
brainmets.freq.num$Normal.brain <- freq.normal
brainmets.freq.num$Gain.brain <- freq.gain

# Fisher tests
loss.pv <- gain.pv <- NULL
for(i in 1:nrow(brainmets.freq.num)) {
  if(i%%100==0) {
    cat(i,"\n")
  }
  # Loss
  mat <- matrix(c(brainmets.freq.num$Loss.primary[i],brainmets.freq.num$Loss.brain[i],sum(brainmets.freq.num$Normal.primary[i],brainmets.freq.num$Gain.primary[i]),sum(brainmets.freq.num$Normal.brain[i],brainmets.freq.num$Gain.brain[i])),byrow = T,nrow = 2,ncol = 2, dimnames = list(c("Loss","Not.Loss"),c("primary","brain")))
  loss.pv <- c(loss.pv,fisher.test(mat, alternative = "two.sided")$p.value)
  # Gain
  mat <- matrix(c(brainmets.freq.num$Gain.primary[i],brainmets.freq.num$Gain.brain[i],sum(brainmets.freq.num$Normal.primary[i],brainmets.freq.num$Loss.primary[i]),sum(brainmets.freq.num$Normal.brain[i],brainmets.freq.num$Loss.brain[i])),byrow = T,nrow = 2,ncol = 2, dimnames = list(c("Gain","Not.Gain"),c("primary","brain")))
  gain.pv <- c(gain.pv,fisher.test(mat, alternative = "two.sided")$p.value)
}
loss.pv.fdr <- p.adjust(loss.pv,method = "fdr")
gain.pv.fdr <- p.adjust(gain.pv,method = "fdr")
summary(loss.pv.fdr)
summary(gain.pv.fdr)
brainmets.freq.num$Loss.pv <- loss.pv
brainmets.freq.num$Loss.pv.fdr <- loss.pv.fdr
brainmets.freq.num$Gain.pv <- gain.pv
brainmets.freq.num$Gain.pv.fdr <- gain.pv.fdr

# Adding frequencies for each group
# primary
freq.loss <- freq.normal <- freq.gain <- NULL
for(i in 1:nrow(brainmets.freq.num)) {
  if(i%%1000==0) {
    cat(i,"\n")
  }
  freq.loss <- c(freq.loss,round(c(brainmets.freq.num$Loss.primary[i]/11)*100))
  freq.normal <- c(freq.normal,round(c(brainmets.freq.num$Normal.primary[i]/11)*100))
  freq.gain <- c(freq.gain,round(c(brainmets.freq.num$Gain.primary[i]/11)*100))
}
brainmets.freq.num$Loss.primary.freq <- freq.loss
brainmets.freq.num$Normal.primary.freq <- freq.normal
brainmets.freq.num$Gain.primary.freq <- freq.gain
# Study Entry
freq.loss <- freq.normal <- freq.gain <- NULL
for(i in 1:nrow(brainmets.freq.num)) {
  if(i%%1000==0) {
    cat(i,"\n")
  }
  freq.loss <- c(freq.loss,round(c(brainmets.freq.num$Loss.brain[i]/13)*100))
  freq.normal <- c(freq.normal,round(c(brainmets.freq.num$Normal.brain[i]/13)*100))
  freq.gain <- c(freq.gain,round(c(brainmets.freq.num$Gain.brain[i]/13)*100))
}
brainmets.freq.num$Loss.brain.freq <- freq.loss
brainmets.freq.num$Normal.brain.freq <- freq.normal
brainmets.freq.num$Gain.brain.freq <- freq.gain

# Frequency plots without Fisher test significance
GroupName <- "BrainMets"
freq.tab <- brainmets.freq.num
freq.tab$Loss.primary.freq <- freq.tab$Loss.primary.freq*-1
freq.tab$Gain.primary.freq <- freq.tab$Gain.primary.freq*-1
cairo_pdf(paste("FrequencyPlots_All_CNtypes_",GroupName,"_woFisherSignificance_v2_COSMIC.pdf",sep=""), width = 14, height = 10, onefile=TRUE)
par(mfrow=c(2,1))

# GAIN
plot(x=xlimit,y=c(-100,100),type="n", xaxt="n", yaxt="n",xlab = "", ylab="Percentage",main = "Gain",cex = 1.5, cex.lab=1.5, cex.main = 1.5)
axis(side=2, at=seq(-100,100,by=25), cex.axis=1.5,las=1,labels=c(100,75,50,25,0,25,50,75,100))
axis(side=1, at=meanind, cex.axis=1,las=1,labels=chr,lwd.ticks=1)
for(i in 1:nrow(freq.tab)) {
  if(freq.tab$chromosome[i]=="1") {
    colo <- "cornflowerblue"
    rect(freq.tab$start[i],0,freq.tab$end[i],freq.tab$Gain.primary.freq[i],col=colo,border=NA)
    rect(freq.tab$start[i],0,freq.tab$end[i],freq.tab$Gain.brain.freq[i],col=colo,border=NA)
  }
  else {
    colo <- "cornflowerblue"
    start <- freq.tab$start[i]+ref$LengthAdded[which(ref$Chromosome==freq.tab$chromosome[i])-1]
    end <- freq.tab$end[i]+ref$LengthAdded[which(ref$Chromosome==freq.tab$chromosome[i])-1]
    rect(start,0,end,freq.tab$Gain.primary.freq[i],col=colo,border=NA)
    rect(start,0,end,freq.tab$Gain.brain.freq[i],col=colo,border=NA)
  }
}
# abline(v=allind,lty="dotted",lwd=0.3,a=1,b=2)
abline(h=0,lty="solid",lwd=1)
abline(h=c(-100,-75,-50,-25,0,25,50,75,100),lty="dotted",lwd=0.3)
abline(v=ref$End,lty="dotted",lwd=0.3)
legend("topleft",legend = "brain",cex = 1.5)
legend("bottomleft",legend = "primary",cex = 1.5)

# LOSS
plot(x=xlimit,y=c(-100,100),type="n", xaxt="n", yaxt="n",xlab = "", ylab="Percentage",main = "Loss",cex = 1.5, cex.lab=1.5, cex.main = 1.5)
axis(side=2, at=seq(-100,100,by=25), cex.axis=1.5,las=1,labels=c(100,75,50,25,0,25,50,75,100))
axis(side=1, at=meanind, cex.axis=1,las=1,labels=chr,lwd.ticks=1)
for(i in 1:nrow(freq.tab)) {
  if(freq.tab$chromosome[i]=="1") {
    colo <- "firebrick3"
    rect(freq.tab$start[i],0,freq.tab$end[i],freq.tab$Loss.primary.freq[i],col=colo,border=NA)
    rect(freq.tab$start[i],0,freq.tab$end[i],freq.tab$Loss.brain.freq[i],col=colo,border=NA)
  }
  else {
    colo <- "firebrick"
    start <- freq.tab$start[i]+ref$LengthAdded[which(ref$Chromosome==freq.tab$chromosome[i])-1]
    end <- freq.tab$end[i]+ref$LengthAdded[which(ref$Chromosome==freq.tab$chromosome[i])-1]
    rect(start,0,end,freq.tab$Loss.primary.freq[i],col=colo,border=NA)
    rect(start,0,end,freq.tab$Loss.brain.freq[i],col=colo,border=NA)
  }
}
abline(h=0,lty="solid",lwd=1)
abline(h=c(-100,-75,-50,-25,0,25,50,75,100),lty="dotted",lwd=0.3)
abline(v=ref$End,lty="dotted",lwd=0.3)
legend("topleft",legend = "Brain",cex = 1.5)
legend("bottomleft",legend = "primary",cex = 1.5)
dev.off()

######################################
########segment to gene CN############
######################################
setwd("C:/Users/hecro/OneDrive - Imperial College London/Brain Metastases in Ovarian Cancer/BrainMets_CNS2/data/100kb/")
load(file="tab.lis.ace.RData")
# t(as.data.frame(lapply(tab.lis,dim)))

library(biomaRt)
library(GenomicRanges,quietly=T)
library(reshape2)

# Loading latest gene annotations for GRCh37 (hg19) (already downloaded from ENSEMBL)
# grch37 <- read.table("~/Box/GeneralData/Genes/hg19/hg19_all_genes_28-09-2021_HM.txt",header = T,sep = "\t",stringsAsFactors = F)
grch37 <- read.table("hg19_all_genes_28-09-2021_HM.txt",header = T,sep = "\t",stringsAsFactors = F)
grch37 <- grch37[,c(1,6,4,5,9,8,7)]
colnames(grch37) <- c("ENSG","CHROM","START","END","HGNC","BAND","STRAND")

grch37.gr <- GRanges(seqnames=grch37$CHROM,ranges = IRanges(grch37$START, grch37$END),ENSG=grch37$ENSG,HGNC=grch37$HGNC,BAND=grch37$BAND,STRAND=grch37$STRAND)

cn.df <- grch37
abs.cn <- grch37

load("final_table_everything_HM.RData")

for(i in 1:length(tab.lis)) {
  cat("Working with ",names(tab.lis)[i], " | ",i,"/",length(tab.lis),"\n")
  file <- tab.lis[[i]]
  colnames(file)[ncol(file)] <- "aCN"
  pl <- final$ACE_Ploidy[which(final$IGF==names(tab.lis)[i])]
  # Categorising CN as per COSMIC guidelines
  # https://cancer.sanger.ac.uk/cosmic/help/cnv/overview
  # Ploidy <= 2.7
  # Gain:
  # average genome ploidy <= 2.7 AND total copy number >= 5
  # Loss:
  #   average genome ploidy <= 2.7 AND total copy number = 0
  # Ploidy > 2.7
  # Gain
  # total copy number >= 9
  # Loss
  # total copy number < ( average genome ploidy - 2.7 )
  if(pl<=2.7) {
    file$CNtype <- "Normal"
    file$CNtype[which(file$aCN<1)]  <- "Loss"
    file$CNtype[which(file$aCN>=5)]  <- "Gain"
  } else if(pl>2.7) {
    magicnumber <- pl-2.7
    file$CNtype <- "Normal"
    file$CNtype[which(file$aCN<magicnumber)]  <- "Loss"
    file$CNtype[which(file$aCN>=9)]  <- "Gain"
  }
  
  fileG <- GRanges(seqnames = file$chromosome,ranges = IRanges(file$start, file$end), Cn=file$aCN,CNtype=file$CNtype)
  temp <- findOverlaps(grch37.gr, fileG, type="within")
  cn.df$Cn <- NA
  cn.df$Cn[queryHits(temp)] <- file$CNtype[subjectHits(temp)]
  colnames(cn.df)[ncol(cn.df)] <- names(tab.lis)[i]
  abs.cn$Cn <- NA
  abs.cn$Cn[queryHits(temp)] <- file$aCN[subjectHits(temp)]
  colnames(abs.cn)[ncol(abs.cn)] <- names(tab.lis)[i]
}

write.csv(cn.df,file="Genes_CNtypes_All_Samples_HM_CosmicDefinition.csv",row.names = F,quote = F)
write.csv(abs.cn,file="Genes_aCN_All_Samples_HM_CosmicDefinition.csv",row.names = F,quote = F)

# Selected genes from AmpliSeq panel
panel <- read.csv("IAA21043_167_coverage_details_HM.csv",header = T,stringsAsFactors = F)
genes <- unique(panel$Gene_Symbol)
# Adding more genes, some were in ZC's panel, some from literature (GG)
genes <- c(genes,"MYC","MECOM","CCNE1","TERT","KRAS","CDKN2A","CDKN2B","CCND2","CCND3")
all(genes %in% cn.df$HGNC) # TRUE!

# Top genes with most copy number abberations across samples
tempdf <- cn.df
tempdf <- tempdf[-which(duplicated(tempdf$ENSG)),]
rownames(tempdf) <- tempdf$ENSG
tempdf <- tempdf[,8:ncol(tempdf)]
cnchange <- NULL
for(i in 1:nrow(tempdf)) {
  if(i%%1000==0) {
    cat(i,"\n")
  }
  cnchange <- c(cnchange,length(which(tempdf[i,]=="Gain" | tempdf[i,]=="Loss")))
}
# Top 20 highly abberated genes
cn.df.selected <- cn.df[match(rownames(tempdf)[1:20],cn.df$ENSG),]
write.csv(cn.df.selected, file = "top_20_highly_abberated_genes", row.names = FALSE)

cn.df.selected <- cn.df[which(cn.df$HGNC %in% genes),]
write.csv(cn.df.selected, file = "abberated_genes_in_octopus", row.names = FALSE)

# Preparing input for ONCOPRINT from: https://www.cbioportal.org/oncoprinter
# load("final_table_everything_HC.RData")
primary <- final$IGF[which(final$SampleType=="primary")]
brain <- final$IGF[which(final$SampleType=="brain")]

cn.df.selected.primary <- cn.df.selected[,c(5,which(colnames(cn.df.selected) %in% primary))]
cn.df.selected.primary <- melt(cn.df.selected.primary,"HGNC")
cn.df.selected.primary$whatever <- "CNA"
cn.df.selected.primary <- cn.df.selected.primary[,c(2,1,3,4)]
#cn.df.selected.primary <- cn.df.selected.primary[-which(is.na(cn.df.selected.primary$value)),]
cn.df.selected.primary <- cn.df.selected.primary[-which(cn.df.selected.primary$value=="Normal"),]
cn.df.selected.primary$value[which(cn.df.selected.primary$value=="Gain")] <- "AMP"
cn.df.selected.primary$value[which(cn.df.selected.primary$value=="Loss")] <- "HOMDEL"
# cn.df.selected.primary.2tier <- cn.df.selected.primary[which(cn.df.selected.primary$value=="AMP" | cn.df.selected.primary$value=="HOMDEL"),]
write.table(cn.df.selected.primary,"CN_primary_for_ONCOPRINT_HM_v2.txt",row.names = F,col.names = F,sep = "\t",quote = F)
# write.table(cn.df.selected.primary.2tier,"CN_primary_for_ONCOPRINT_2T_HM_v2.txt",row.names = F,col.names = F,sep = "\t",quote = F)

cn.df.selected.brain <- cn.df.selected[,c(5,which(colnames(cn.df.selected) %in% brain))]
cn.df.selected.brain <- melt(cn.df.selected.brain,"HGNC")
cn.df.selected.brain$whatever <- "CNA"
cn.df.selected.brain <- cn.df.selected.brain[,c(2,1,3,4)]
cn.df.selected.brain <- cn.df.selected.brain[-which(is.na(cn.df.selected.brain$value)),]
cn.df.selected.brain <- cn.df.selected.brain[-which(cn.df.selected.brain$value=="Normal"),]
cn.df.selected.brain$value[which(cn.df.selected.brain$value=="Gain")] <- "AMP"
cn.df.selected.brain$value[which(cn.df.selected.brain$value=="Loss")] <- "HOMDEL"
# cn.df.selected.brain.2tier <- cn.df.selected.brain[which(cn.df.selected.brain$value=="AMP" | cn.df.selected.brain$value=="HOMDEL"),]
write.table(cn.df.selected.brain,"CN_brain_for_ONCOPRINT_HM_v2.txt",row.names = F,col.names = F,sep = "\t",quote = F)
# write.table(cn.df.selected.brain.2tier,"CN_brain_for_ONCOPRINT_2T_HM_v2.txt",row.names = F,col.names = F,sep = "\t",quote = F)

genes <- names(sort(table(cn.df.selected.brain$HGNC),decreasing = T))
# MYC,PALB2,PIK3CA,CDKN2A,CDKN2B,CCNE1,IGFLR1,MECOM,AKT2,BRCA2,NF1,AKT3,PTEN,RB1,TERT,BRCA1,CCND1,CCND2,CDK12,KRAS,RAD51B,RAD51D,AKT1,BARD1,FANCM,TP53
# Previously
# MYC,MECOM,PIK3CA,PALB2,CDKN2A,CDKN2B,KRAS,NF1,TERT,AKT1,CCND2,CCNE1,PTEN,RAD51D,AKT3,IGFLR1,AKT2,BRCA1,CCND1,CDK12,BRCA2,CCND3,TP53,BARD1,BRIP1,RB1,RAD51B,RAD51C,RPS6KB1,FANCM 
# Now go to https://www.cbioportal.org/oncoprinter and make oncoprints

# Comparing all genes between primary and study entry
# Fisher tests
primary.for.ft <- cn.df[,c(5,which(colnames(cn.df) %in% primary))]
primary.for.ft <- primary.for.ft[-which(duplicated(cn.df$HGNC)),]
rownames(primary.for.ft) <- primary.for.ft$HGNC
primary.for.ft <- primary.for.ft[,2:ncol(primary.for.ft)]
brain.for.ft <- cn.df[,c(5,which(colnames(cn.df) %in% brain))]
brain.for.ft <- brain.for.ft[-which(duplicated(cn.df$HGNC)),]
rownames(brain.for.ft) <- brain.for.ft$HGNC
brain.for.ft <- brain.for.ft[,2:ncol(brain.for.ft)]
all(rownames(primary.for.ft)==rownames(brain.for.ft)) # TRUE
pvs <- pvs.Loss <- pvs.Gain <- NULL
for(i in 1:nrow(primary.for.ft)) {
  if(i%%1000==0) {
    cat(i,"\n")
  }
  # Both
  var1 <- length(which(primary.for.ft[i,]=="Loss" | primary.for.ft[i,]=="Gain"))
  var2 <- length(which(brain.for.ft[i,]=="Loss" | brain.for.ft[i,]=="Gain"))
  mat <- matrix(c(var1,var2,ncol(primary.for.ft)-var1,ncol(brain.for.ft)-var2),byrow = T,nrow = 2,ncol = 2, dimnames = list(c("Changed","Not.Changed"),c("primary","brain")))
  pvs <- c(pvs,fisher.test(mat, alternative = "two.sided")$p.value)
  # Loss
  var1 <- length(which(primary.for.ft[i,]=="Loss"))
  var2 <- length(which(brain.for.ft[i,]=="Loss"))
  mat <- matrix(c(var1,var2,ncol(primary.for.ft)-var1,ncol(brain.for.ft)-var2),byrow = T,nrow = 2,ncol = 2, dimnames = list(c("Changed","Not.Changed"),c("primary","brain")))
  pvs.Loss <- c(pvs.Loss,fisher.test(mat, alternative = "two.sided")$p.value)
  # Gain
  var1 <- length(which(primary.for.ft[i,]=="Gain"))
  var2 <- length(which(brain.for.ft[i,]=="Gain"))
  mat <- matrix(c(var1,var2,ncol(primary.for.ft)-var1,ncol(brain.for.ft)-var2),byrow = T,nrow = 2,ncol = 2, dimnames = list(c("Changed","Not.Changed"),c("primary","brain")))
  pvs.Gain <- c(pvs.Gain,fisher.test(mat, alternative = "two.sided")$p.value)
}
pvs.fdr <- p.adjust(pvs,method = "fdr")
pvs.Loss.fdr <- p.adjust(pvs.Loss,method = "fdr")
pvs.Gain.fdr <- p.adjust(pvs.Gain,method = "fdr")
summary(pvs)
summary(pvs.fdr)
summary(pvs.Loss)
summary(pvs.Loss.fdr)
summary(pvs.Gain)
summary(pvs.Gain.fdr)

# order.wise <- pvs.homdel
# names(order.wise) <- rownames(primary.for.ft)
# order.wise <- sort(order.wise,decreasing = F)
