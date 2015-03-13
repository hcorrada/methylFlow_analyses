#!/usr/bin/env Rscript
### run ./sam.r par1

### par1 = input directory


#data[1] = "SRR1020537-3006910-15008000"

library("devtools")
path = "~/Projects/methylFlow/src/methylFlow/methylFlowr"
load_all(path)

data <- commandArgs(T)

dirNew <- "/Users/faezeh/Projects/methylFLow"
dirOld <- "/Users/faezeh/Projects/methylFLow_old"
#install_github("rafalib","ririzarr")
library("rafalib")

datadir <- "data/sam"
resultdirNew <- file.path(paste(dirNew,"/exps/sam/",data[1],sep=""))
resultdirOld <- file.path(dirOld,"exps/sam/auto")



#dir <- dirOld
#resultdir <- resultdirOld

dir <- dirNew
resultdir <- resultdirNew


print(resultdir)

#pd <- data.frame(subject=rep(4:6,2), status=c("T", "T", "T","N", "N", "N"))
j <- "3"

#pd$dirname=sprintf("CAP_%s_%d", pd$status, pd$subject)
#pd$dirname = datadir
objs <- lapply(seq(len=1), function(i) {
    curdir <- file.path(dir,datadir)
    read.methylflow.dir(resultdir,has.header=TRUE)
})
#names(objs) <- file.path(dir,datadir)
names(objs) <- data[1]

objs <- lapply(objs, processMethylpats)
filteredObjs <- lapply(objs, mfFilterBy, minComponentCoverage=50, minComponentWidth=86)
names(filteredObjs) <- names(objs)


figdir <- file.path(resultdir,"figs2")

if (!file.exists(figdir)) dir.create(figdir)

pdf(file.path(figdir,"fragment_length.pdf"),height=9,width=9)
#mypar(1,1)

widths <- lapply(filteredObjs, function(obj) width(components(obj)))
ncpgs <- lapply(filteredObjs, function(obj) ncpgs(obj, level="pattern", summary="max"))

#names(ncpgs) <- names(widths) <- gsub("CAP_", "", names(widths))
boxplot(widths, main="reconstructed fragment size")
boxplot(ncpgs, main="number of cpgs in reconstructed fragments")
dev.off()



####
filteredObjs <- lapply(filteredObjs, mfFilterBy, minNumberOfPatterns=1)

npats <- lapply(filteredObjs, npatterns, by.component=TRUE)
compCoverage <- lapply(filteredObjs, counts, level="component", kind="raw")

pdf(file.path(figdir,"patt_by_coverage.pdf"), height=9, width=9)
#mypar(1,1)
  
for (i in seq(along=filteredObjs)) {
  plot(compCoverage[[i]], npats[[i]],main=names(filteredObjs[[i]]),xlab="coverage",ylab="num. patterns")
}
dev.off()



rawCpgs <- lapply(filteredObjs, makeCpgGR, kind="raw")
estimatedCpgs <- lapply(filteredObjs, makeCpgGR, kind="estimated")

#readGRs <- lapply(seq(len=nrow(pd)), function(i) {
#    curdir <- file.path("reads", pd$dirname[i])
#    tab <- read.delim(file.path(curdir, "13.methylation.withsub.tsv"),
#                      header=FALSE,stringsAsFactors=FALSE)
#    names(tab) <- c("id", "start", "width", "cstrand", "methylpat", "subst")
#    tab$end <- tab$start + tab$width - 1
#    tab$chr <- "13"
#    tab$width <- NULL
#    gr <- tab2gr(tab)
#    gr
#})
#
#.parseMethylpats <- function(x) {
#    tmp <- strsplit(x, ",")
#    ncpgs <- ifelse(x=="*", 0, sapply(tmp,length))
#    tmp2 <- lapply(seq(along=x), function(i) {
#        if (ncpgs[i] == 0) return(NULL)
#        strsplit(tmp[[i]], ":")
#    })
#    
#    locs <- lapply(tmp2, function(y) {
#      if (is.null(y)) return(0)
#      as.integer(sapply(y,"[",1))
#    })
#
#    meth <- lapply(tmp2, function(y) {
#        if (is.null(y)) return("")
#        sapply(y,"[",2)
#    })
#    list(ncpgs=ncpgs,locs=locs,meth=meth)
#}
#
#readGRs <- lapply(readGRs, function(x) {
#    tmp <- .parseMethylpats(x$methylpat)
#    x$ncpgs <- tmp$ncpgs
#    x$locs <- tmp$locs
#    x$meth <- tmp$meth
#    x
#})
#
#readCpgs <- lapply(readGRs, function(gr) {
#    keep <- gr$ncpgs > 0
#    gr <- gr[keep,]
#    locs <- rep(start(gr), gr$ncpgs) + unlist(gr$locs)
#
#    Cov <- rep(1, length(locs))
#    Meth <- 1*(unlist(gr$meth) == "M")
#
#    Cov <- tapply(Cov, locs, sum)
#    Meth <- tapply(Meth, locs, sum)
#    Loc <- tapply(locs, locs, function(x) x[1])
#    chr <- tapply(rep(as.character(seqnames(gr)), gr$ncpgs), locs, function(x) x[1])
#    newGR <- GRanges(chr, IRanges(start=Loc, width=1), Cov=Cov, Meth=Meth)
#    newGR$Beta <- Meth / Cov
#    names(newGR) <- NULL
#    newGR
#})
#
#readCpgs <- lapply(readCpgs, renameSeqlevels, value=structure("chr13", names="13"))
#
betaGR <- lapply(seq(along=rawCpgs), function(i) {
    rawGR <- rawCpgs[[i]]
    estGR <- estimatedCpgs[[i]]
#    readGR <- readCpgs[[i]]

    olaps <- findOverlaps(rawGR, estGR)
    gr <- rawGR[queryHits(olaps),]
    out <- gr
    mcols(out) <- NULL
    
    out$rawBeta <- gr$Beta
    out$rawCov <- gr$Cov
    
    gr <- estGR[subjectHits(olaps),]
    out$estimatedBeta <- gr$Beta

#    olaps <- findOverlaps(readGR, out)
 #   out$readBeta <- readGR$Beta[queryHits(olaps)]
    out
})

pdf(file.path(figdir,"meth_percentage.pdf"), height=9, width=9)
#mypar(1,1)

for (i in seq(along=betaGR)) {
    gr <- betaGR[[i]]
    plot(gr$rawBeta,gr$estimatedBeta,
         bty='l',
         main=names(filteredObjs)[i],
         ylab="pattern methyl Percentage",
         xlab="region methyl Precentage",
         #cex = 0.2,
         cex.lab=.9)
}

#for (i in seq(along=betaGR)) {
#    gr <- betaGR[[i]]
#    plot(gr$rawBeta,gr$readBeta,
#         bty='l',
#         main=names(filteredObjs)[i],
#         ylab="read methyl Percentage",
#         xlab="region methyl Precentage",
#         cex = 0.2,
#         cex.lab=.9)
#}
#
dev.off()

