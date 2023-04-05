.makeSegments <- function(data,chrdata) {
  previous    <- 2000
  chrpr       <- -100
  values      <- c()
  start       <- c()
  end         <- c()
  el          <- length(data)
  data <- c(data,-10000) #add value to allow data[i+1]
  for (i in 1:el) {
    if ((data[i] != previous & previous != data[i+1]) | chrdata[i] != chrpr) { #bug repaired 12/06/09
      start   <- c(start, i)
      last    <- i - 1
      if (last > 0) end <- c(end, last)
      values  <- c(values, data[i])
    }
    previous    <- data[i]
    chrpr <- chrdata[i]
  }
  end     <- c(end, el)
  result  <- cbind(values, start, end)
  result
}

.getChromosomeLengths <- function(build) {
  build <- as.integer(gsub('[^0-9]', '', build))
  if (build == 34 || build == 16) {
    chromosome.lengths <- c(246127941, 243615958, 199344050, 191731959, 181034922, 170914576, 158545518, 146308819, 136372045, 135037215, 134482954, 132078379, 113042980, 105311216, 100256656, 90041932, 81860266, 76115139, 63811651, 63741868, 46976097, 49396972, 153692391, 50286555)
  } else if (build == 35 || build == 17) {
    chromosome.lengths <- c(245522847, 243018229, 199505740, 191411218, 180857866, 170975699, 158628139, 146274826, 138429268, 135413628, 134452384, 132449811, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49554710, 154824264, 57701691)
  } else if (build == 36 || build == 18) {
    chromosome.lengths <- c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else {
    chromosome.lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  }
  names(chromosome.lengths) <- 1:24
  chromosome.lengths
}

#########################################################################/**
# @RdocFunction plot
#
# @alias plot,QDNAseqSignals,missing-method
#
# @title "Plot copy number profile"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{x}{A @see "QDNAseqReadCounts" or @see "QDNAseqCopyNumbers" object.}
#     \item{y}{missing}
#     \item{...}{...}
#%     \item{verbose}{If @TRUE, verbose messages are produced.}
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# copyNumbers <- correctBins(readCountsFiltered)
# plot(copyNumbers)
# }
#
# @author "IS"
#
# @keyword hplot
#*/#########################################################################
plot_fix <- function (x, y, main=NULL, includeReadCounts=TRUE,
                    logTransform=TRUE, scale=TRUE, sdFUN="sdDiff",
                    delcol=getOption("QDNAseq::delcol", "darkred"),
                    losscol=getOption("QDNAseq::losscol", "red"),
                    gaincol=getOption("QDNAseq::gaincol", "blue"),
                    ampcol=getOption("QDNAseq::ampcol", "darkblue"),
                    pointcol=getOption("QDNAseq::pointcol", "black"),
                    segcol=getOption("QDNAseq::segcol", "chocolate"),
                    misscol=getOption("QDNAseq::misscol", NA),
                    pointpch=getOption("QDNAseq::pointpch", 1L),
                    pointcex=getOption("QDNAseq::pointcex", 0.1),
                    xlab=NULL, ylab=NULL, ylim=NULL, xaxt="s", yaxp=NULL,
                    showDataPoints=TRUE, showSD=TRUE, doSegments=TRUE, doCalls=TRUE, ...,
                    verbose=getOption("QDNAseq::verbose", TRUE)) {
            
            
            if (inherits(x, c("QDNAseqCopyNumbers", "QDNAseqReadCounts"))) {
              condition <- binsToUse(x)
            } else {
              condition <- rep(TRUE, times=nrow(x))
            }
            baseLine <- NA_real_
            doCalls <- "calls" %in% assayDataElementNames(x) & doCalls
            doSegments <- "segmented" %in% assayDataElementNames(x) & doSegments
            if (doCalls) {
              if (is.null(ylim))
                if (logTransform) {
                  ylim <- c(-5, 5)
                } else {
                  ylim <- c(-2, 4)
                }
            }
            if ("copynumber" %in% assayDataElementNames(x)) {
              copynumber <- assayDataElement(x, "copynumber")[condition, , drop=FALSE]
              if (is.null(ylab))
                ylab <- ifelse(logTransform, expression(log[2]~ratio), "ratio")
              if (is.null(ylim))
                if (logTransform) {
                  ylim <- c(-3, 5)
                } else {
                  ylim <- c(0, 4)
                }
              if (is.null(yaxp))
                yaxp <- c(ylim[1], ylim[2], ylim[2]-ylim[1])
              baseLine <- ifelse(logTransform, 0, 1)
            } else {
              copynumber <- assayDataElement(x, "counts")[condition, , drop=FALSE]
              if (is.null(ylab))
                ylab <- ifelse(logTransform, expression(log[2]~read~count),
                               "read count")
              if (is.null(ylim))
                if (logTransform) {
                  ylim <- c(0, max(log2adhoc(copynumber)))
                } else {
                  ylim <- range(copynumber)
                }
            }
            if (is.null(main))
              main <- sampleNames(x)
            if (includeReadCounts && "total.reads" %in% names(pData(x)))
              main <- paste(main, " (",
                            format(x$total.reads, trim=TRUE, big.mark=","), " reads)", sep="")
            if (length(ylab) == 1)
              ylab <- rep(ylab, times=ncol(x))
            all.chrom <- chromosomes(x)
            if (is.integer(all.chrom)) # when x is a cghRaw, cghSeg, or cghCall object
              all.chrom <- as.character(all.chrom)
            chrom <- all.chrom[condition]
            uni.chrom <- unique(chrom)
            chrom.num <- as.integer(factor(chrom, levels=uni.chrom, ordered=TRUE))
            uni.chrom.num <- unique(chrom.num)
            if (!scale) {
              pos <- pos2 <- 1:sum(condition)
              chrom.ends <- aggregate(pos,
                                      by=list(chromosome=chrom), FUN=max)$x
            } else {
              if (inherits(x, c("cghRaw", "cghSeg", "cghCall"))) {
                chrom.lengths <- .getChromosomeLengths("GRCh37")
              } else {
                all.chrom.lengths <- aggregate(bpend(x),
                                               by=list(chromosome=all.chrom), FUN=max)
                chrom.lengths <- all.chrom.lengths$x
                names(chrom.lengths) <- all.chrom.lengths$chromosome
              }
              pos <- as.numeric(bpstart(x)[condition])
              pos2 <- as.numeric(bpend(x)[condition])
              chrom.lengths <- chrom.lengths[uni.chrom]
              chrom.ends <- integer()
              cumul <- 0
              for (i in seq_along(uni.chrom)) {
                pos[chrom.num > uni.chrom.num[i]] <-
                  pos[chrom.num > uni.chrom.num[i]] +
                  chrom.lengths[uni.chrom[i]]
                pos2[chrom.num > uni.chrom.num[i]] <-
                  pos2[chrom.num > uni.chrom.num[i]] +
                  chrom.lengths[uni.chrom[i]]
                cumul <- cumul + chrom.lengths[uni.chrom[i]]
                chrom.ends <- c(chrom.ends, cumul)
              }
              names(chrom.ends) <- names(chrom.lengths)
            }
            if (length(uni.chrom) == 1) {
              xax <- pretty(c(0, chrom.lengths[uni.chrom]))
              xaxlab <- xax / 1e6L
              if (is.null(xlab))
                xlab <- paste0("chromosome ", uni.chrom, ", Mbp")
            } else {
              xax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
              xaxlab <- uni.chrom
              if (is.null(xlab))
                xlab <- "chromosome"
            }
            if (inherits(x, c("cghRaw", "cghSeg", "cghCall")))
              copynumber <- log2adhoc(copynumber, inv=TRUE)
            #if (is.character(sdFUN) && sdFUN == "sdDiff") {
              symbol <- quote(hat(sigma)[Delta^"*"])
            #} else if (is.character(sdFUN) && length(grep("Diff", sdFUN)) == 1) {
            #  symbol <- quote(hat(sigma)[Delta])
            #} else {
            #  symbol <- quote(hat(sigma))
            #}
            sdFUN <- match.fun(sdFUN)
            noise <- apply(copynumber, MARGIN=2L, FUN=sdFUN, na.rm=TRUE)
            if (logTransform)
              copynumber <- log2adhoc(copynumber)
            for (i in seq_len(ncol(x))) {
              # vmsg("Plotting sample ", main[i], " (", i, " of ", ncol(x), ") ...",
              #      appendLF=FALSE)
              cn <- copynumber[, i]
              if (doSegments) {
                segmented <- assayDataElement(x, "segmented")[condition, i]
                if (inherits(x, c("cghRaw", "cghSeg", "cghCall")))
                  segmented <- log2adhoc(segmented, inv=TRUE)
                if (logTransform)
                  segmented <- log2adhoc(segmented)
                segment <- .makeSegments(segmented, chrom)
              }
              if (doCalls) {
                losses <- probloss(x)[condition, i]
                gains <- probgain(x)[condition, i]
                if (!is.null(probdloss(x)))
                  losses <- losses + probdloss(x)[condition, i]
                if (!is.null(probamp(x)))
                  gains <- gains + probamp(x)[condition, i]
                par(mar=c(5, 4, 4, 4) + 0.2)
                plot(NA, main=main[i], xlab=NA, ylab=NA, las=1,
                     xlim=c(0, max(chrom.ends)), ylim=ylim, xaxs="i", xaxt="n",
                     yaxp=c(ylim[1], ylim[2], ylim[2]-ylim[1]), yaxs="i")
                lim <- par("usr")
                lim[3:4] <- c(0, 1)
                par(usr=lim)
                dticks <- seq(0, 1, by=0.2)
                axis(4, at=dticks, labels=NA, srt=270, las=1, cex.axis=1,
                     cex.lab=1, tck=-0.015)
                axis(4, at=dticks, labels=dticks, srt=270, las=1, cex.axis=1,
                     cex.lab=1, line=-0.4, lwd=0)
                mtext("probability", side=4, line=2, cex=par("cex"))
                if (!is.na(misscol)) {
                  rect(0, -1, max(pos2), 1, col=misscol, border=NA)
                  rect(pos, -1, pos2, 1, col="white", border=NA)
                }
                rect(pos[segment[,2]], 0, pos2[segment[,3]], losses[segment[,2]],
                     col=losscol, border=losscol)
                if (!is.null(probdloss(x)))
                  rect(pos[segment[,2]], 0, pos2[segment[,3]],
                       probdloss(x)[condition, i][segment[,2]],
                       col=delcol, border=delcol)
                rect(pos[segment[,2]], 1, pos2[segment[,3]], 1-gains[segment[,2]],
                     col=gaincol, border=gaincol)
                if (!is.null(probamp(x)))
                  rect(pos[segment[,2]], 1, pos2[segment[,3]],
                       1-probamp(x)[condition, i][segment[,2]],
                       col=ampcol, border=ampcol)
                axis(3, at=pos[which(probamp(x)[condition,i] >= 0.5)],
                     labels=FALSE, col=ampcol, col.axis="black", srt=270, las=1,
                     cex.axis=1, cex.lab=1)
                axis(1, at=pos[which(probdloss(x)[condition,i] >= 0.5)],
                     labels=FALSE, col=delcol, col.axis="black", srt=270, las=1,
                     cex.axis=1, cex.lab=1)
                box()
                lim[3:4] <- ylim
                par(usr=lim)
                points(pos, cn, cex=pointcex, col=pointcol, pch=pointpch)
              } else {
                plot(pos, cn, cex=pointcex, col=pointcol, main=main[i],
                     xlab=NA, ylab=NA, ylim=ylim, xaxt="n", xaxs="i", yaxs="i",
                     yaxp=yaxp, tck=-0.015, las=1, pch=pointpch)
              }
              mtext(text=xlab, side=1, line=2, cex=par("cex"))
              mtext(text=ylab[i], side=2, line=2, cex=par("cex"))
              abline(h=baseLine)
              abline(v=chrom.ends[-length(chrom.ends)], lty="dashed")
              if (!is.na(xaxt) && xaxt != "n") {
                axis(side=1, at=xax, labels=NA, cex=.2, lwd=.5, las=1,
                     cex.axis=1, cex.lab=1, tck=-0.015)
                axis(side=1, at=xax, labels=xaxlab, cex=.2, lwd=0, las=1,
                     cex.axis=1, cex.lab=1, tck=-0.015, line=-0.4)
              }
              if (doSegments) {
                for (jjj in seq_len(nrow(segment))) {
                  segments(pos[segment[jjj,2]], segment[jjj,1],
                           pos[segment[jjj,3]], segment[jjj,1], col=segcol, lwd=3)
                }
              }
              par(xpd=TRUE)
              amps <- cn
              amps[amps <= ylim[2]] <- NA_real_
              amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
              dels <- cn
              dels[dels >= ylim[1]] <- NA_real_
              dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
              points(pos, amps, pch=24, col=pointcol, bg=pointcol, cex=0.5)
              points(pos, dels, pch=25, col=pointcol, bg=pointcol, cex=0.5)
              if (doSegments) {
                amps <- assayDataElement(x, "segmented")[condition, i]
                if (logTransform)
                  amps <- log2adhoc(amps)
                amps[amps <= ylim[2]] <- NA_real_
                amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
                dels <- assayDataElement(x, "segmented")[condition, i]
                if (logTransform)
                  dels <- log2adhoc(dels)
                dels[dels >= ylim[1]] <- NA_real_
                dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
                points(pos, amps, pch=24, col=segcol, bg=segcol, cex=0.5)
                points(pos, dels, pch=25, col=segcol, bg=segcol, cex=0.5)
              }
              par(xpd=FALSE)
              ### estimate for standard deviation
              if (showSD) {
                if (!is.na(x$expected.variance[i])) {
                  sdexp <- substitute(paste(E~sigma==e, ", ", symbol==sd),
                                      list(e=sprintf("%.3g", sqrt(x$expected.variance[i])),
                                           symbol=symbol, sd=sprintf("%.3g", noise[i])))
                } else {
                  sdexp <- substitute(symbol==sd,
                                      list(symbol=symbol, sd=sprintf("%.3g", noise[i])))
                }
                mtext(sdexp, side=3, line=0, adj=1, cex=par("cex"))
              }
              ### number of data points
              if (showDataPoints) {
                str <- paste(round(sum(condition) / 1000), "k x ", sep="")
                probe <- median(bpend(x)-bpstart(x)+1)
                if (probe < 1000) {
                  str <- paste(str, probe, " bp", sep="")
                } else {
                  str <- paste(str, round(probe / 1000), " kbp", sep="")
                }
                if (doSegments)
                  str <- paste(str, ", ", nrow(segment), " segments", sep="")
                mtext(str, side=3, line=0, adj=0, cex=par("cex"))
              }
              return(nrow(segment))
            }
}

sdDiffTrim <- function(x, ..., trim=0.001, scale=TRUE) {
  if (scale)
    x <- x / mean(x, na.rm=TRUE)
  matrixStats::sdDiff(x, ..., trim=trim)
}

expectedVariance <- function(object) {
  expected.variance <- sum(binsToUse(object)) / object$used.reads
  if ("paired.ends" %in% colnames(pData(object))) {
    multiplier <- ifelse(object$paired.ends, 2, 1)
    expected.variance <- expected.variance * multiplier
  }
  expected.variance
}

log2offset <- function(offset=.Machine$double.xmin) {
  # Get offset?
  if (missing(offset)) {
    offset <- getOption("QDNAseq::log2offset", .Machine$double.xmin)
    offset <- as.double(offset);
    stopifnot(is.finite(offset));
    return(offset);
  }
  
  # Reset offset?
  if (is.null(offset)) offset <- eval(formals(log2offset)$offset);
  
  # Set offset
  stopifnot(length(offset) == 1L);
  offset <- as.double(offset);
  stopifnot(is.finite(offset));
  options("QDNAseq::log2offset"=offset);
  
  offset;
}

log2adhoc <- function(x, offset=log2offset(), inv=FALSE) {
  if (!inv) {
    x[x < 0] <- 0
    x <- x + offset
    log2(x)
  } else {
    x <- 2^x
    x - offset
  }
}
#### do not touch!!!!


library(matrixStats)
library(QDNAseq)
library(purrr)


get_segments <- function(x) {
  condition <- binsToUse(x)
  if (inherits(x, c("cghRaw", "cghSeg", "cghCall")))
    segmented <- log2adhoc(segmented, inv=TRUE)
  if (logTransform)
    segmented <- log2adhoc(segmented)
  segment <- .makeSegments(segmented, chrom)
  
  copynumber <- assayDataElement(x, "copynumber")[condition, , drop=FALSE]
  doSegments=T
  logTransform=T
  all.chrom <- chromosomes(x)
  if (is.integer(all.chrom)) # when x is a cghRaw, cghSeg, or cghCall object
    all.chrom <- as.character(all.chrom)
  chrom <- all.chrom[condition]
  
  1:length(colnames(copynumber)) %>%
    map(function(i) {
      cn <- copynumber[, i]
      segmented <- assayDataElement(x, "segmented")[condition, i]
      if (inherits(x, c("cghRaw", "cghSeg", "cghCall")))
        segmented <- log2adhoc(segmented, inv=TRUE)
      if (logTransform)
        segmented <- log2adhoc(segmented)
      segment <- .makeSegments(segmented, chrom)
      tibble::tibble(n_segments = nrow(segment))
    }) %>%
    setNames(colnames(copynumber)) %>%
    dplyr::bind_rows(.id = 'sample')
}
