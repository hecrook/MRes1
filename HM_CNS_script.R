#R code script 2
#######################COSMIC SIGNATURES#####################################
cat("Getting all breakpoints..\n")
Breakpoints <- vector("list", 22)
chr <- c(1:22)
for (i in 1:length(tab.lis)) {
  cat("Working with ",names(tab.lis)[i]," | ",i,"/",length(tab.lis),"\n")
  file <- tab.lis[[i]]
  file$start <- as.numeric(file$start)
  file$end <- as.numeric(file$end)
  for (j in 1:22) {
    vec <- 0
    ind <- which(file$chromosome==chr[j])
    vec<- c(vec,file$start[ind])
    Breakpoints[[j]] <- c(Breakpoints[[j]],vec)
  }
}

#Sorting all breakpoints per chromosome after removing redundant breakpoints
cat("Sorting breakpoints per chromosome and removing duplicates..\n")
for(j in 1:length(Breakpoints)) {
  Breakpoints[[j]] <- sort(unique(Breakpoints[[j]]))
}

#Creating a data fram and turning all breakpoints into strat and end of segments i.e. creating consensus segments
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

#Mapping CN from all samples to consensus segments 
cat("Mapping CN from all TAPS CN files to consensus segments..\n")
# consensusG <- GRanges(seqnames = gsub("chr", "", df$chromosome, fixed=TRUE),ranges = IRanges(df$start, df$end))
consensusG <- GRanges(seqnames = df$chromosome,ranges = IRanges(df$start, df$end))
count <- 0
df2 <- df[,1:3]
for(i in 1:length(tab.lis.2)) {
  cat("Working with ",names(tab.lis.2)[i], " | ",i,"/",length(tab.lis.2),"\n")
  file <- tab.lis.2[[i]]
  file$start <- as.numeric(file$start)
  file$end <- as.numeric(file$end)
  colnames(file)[ncol(file)] <- "aCN"
  file$aCN <- as.character(file$aCN)
  fileG <- GRanges(seqnames = file$chromosome,ranges = IRanges(file$start, file$end), Cn=file$aCN)
  temp <- findOverlaps(consensusG, fileG, type="within")
  
  df2$Cn[queryHits(temp)] <- file$aCN[subjectHits(temp)]
  colnames(df2)[ncol(df2)] <- names(tab.lis.2)[i]
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

OHA.consensus <- df2

#### CALCULATING FREQUENCIES PER GROUP
#Why do we minus 4 from the number of columns??
octopus.cat <- as.data.frame(matrix(nrow = nrow(OHA.consensus),ncol = ncol(OHA.consensus)-4))
colnames(octopus.cat) <- colnames(OHA.consensus)[5:ncol(OHA.consensus)]
# Rounding off to 2 digits
temp <- OHA.consensus[,5:ncol(OHA.consensus)]
for(i in 1:ncol(temp)) {
  temp[,i] <- round(as.numeric(temp[,i]),2)
}

# Upscaling and downscaling values to closest whole value if near than 10% 
for(i in 1:ncol(temp)) {
  temp[which(temp[,i]>=1.9 & temp[,i]<2),i] <- 2
  temp[which(temp[,i]>=2.9 & temp[,i]<3),i] <- 3
  temp[which(temp[,i]>=3.9 & temp[,i]<4),i] <- 4
}
