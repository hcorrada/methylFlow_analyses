library("devtools")
path = "../../src/methylFlow/methylFlowr"
load_all(path)

#install_github("rafalib","ririzarr")
library("rafalib")

datadir <- "cpg_output"
pd <- data.frame(subject=rep(4:6,2), status=c("T", "T", "T","N", "N", "N"))
j <- "13"

pd$dirname=sprintf("CAP_%s_%d", pd$status, pd$subject)
objs <- lapply(seq(len=nrow(pd)), function(i) {
    curdir <- file.path(datadir,pd$dirname[i],j)
    read.methylflow.dir(curdir, pd$dirname[i],has.header=TRUE)
})
names(objs) <- pd$dirname

objs <- lapply(objs, processMethylpats)

filteredObjs <- lapply(objs, mfFilterBy, minComponentCoverage=100, minComponentWidth=86)
names(filteredObjs) <- names(objs)


figdir <- file.path("figs2")
if (!file.exists(figdir)) dir.create(figdir)

pdf(file.path(figdir,"fragment_length.pdf"),height=4,width=4)
mypar(1,1)

widths <- lapply(filteredObjs, function(obj) width(components(obj)))
ncpgs <- lapply(filteredObjs, function(obj) ncpgs(obj, level="pattern", summary="max"))

names(ncpgs) <- names(widths) <- gsub("CAP_", "", names(widths))
boxplot(widths, main="reconstructed fragment size")
boxplot(ncpgs, main="number of cpgs in reconstructed fragments")
dev.off()

####
filteredObjs <- lapply(filteredObjs, mfFilterBy, minNumberOfPatterns=1)

npats <- sapply(filteredObjs, npatterns, by.component=TRUE)
compCoverage <- sapply(filteredObjs, counts, level="component", kind="raw")

pdf(file.path(figdir,"patt_by_coverage.pdf"), height=4, width=6)
mypar(2,3)
  
for (i in seq(along=filteredObjs)) {
    plot(compCoverage[[i]], npats[[i]],main=names(filteredObjs)[i],xlab="coverage",ylab="num. patterns")
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

pdf(file.path(figdir,"meth_percentage.pdf"), height=3, width=3)
mypar(1,1)

for (i in seq(along=betaGR)) {
    gr <- betaGR[[i]]
    plot(gr$rawBeta,gr$estimatedBeta,
         bty='l',
         main=names(filteredObjs)[i],
         ylab="pattern methyl Percentage",
         xlab="region methyl Precentage",
         cex = 0.2,
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


allBeta <- betaGR[[1]]
tmp <- allBeta$estimatedBeta

mcols(allBeta) <- NULL
mcols(allBeta)[[names(filteredObjs)[1]]] <- tmp

for (i in seq(2,length(betaGR))) {
    gr <- betaGR[[i]]
    olaps <- findOverlaps(gr, allBeta)

    sname <- names(filteredObjs)[i]
    mcols(allBeta)[[sname]] <- NA
    mcols(allBeta)[[sname]][subjectHits(olaps)] <- gr$estimatedBeta[queryHits(olaps)]
}

library(bumphunter)
clusters <- clusterMaker(seqnames(allBeta), start(allBeta), maxGap=100)

library(doParallel)
registerDoParallel(cores=2)
betamat <- as.matrix(mcols(allBeta))

bumps <- bumphunter(betamat, model.matrix(~1+status,data=pd),
                    pos=start(allBeta), cluster=clusters, B=0, smooth=FALSE, cutoff=.2)

#plotRegion <- GRanges("chr13", IRanges(start=bumps$tab$start[1], bumps$tab$end[1]))
#plotRegion <- plotRegion * -8

plotRegion <- GRanges("chr13", IRanges(start=38444300,end=38444850))

regionBeta <- subsetByOverlaps(allBeta, plotRegion)


betamat <- as.matrix(mcols(regionBeta))
Index <- splitit(pd$status)
betamns <- sapply(Index, function(ind) rowMeans(betamat[,ind]))

fits <- lapply(1:2, function(i) loess(betamns[,i]~start(regionBeta)))

load_all(path)
pdf(file.path(figdir, "region.pdf"),width=8,height=6)

layoutmat <- matrix(c(8,rep(1,4),9,rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3)),nc=1)
layout(layoutmat)
par(mar=c(0, 3, 0, 1.1))
par(oma=c(0,0,0,0))
    
xlim <- c(start(plotRegion), end(plotRegion))

plot(0, type="n", xlim=xlim, ylim=c(0,1), ylab="",xaxt="n",yaxt="n")
lat <- pretty(xlim,n=5)
axis(side=1, at=lat, cex.axis=.8, mgp=c(1,.3,0))
axis(side=2, at=c(0, 0.5, 1),las=1,cex.axis=.9, mgp=c(1,.4,0))
title(xlab="genomic position (chr13)", ylab="DNA methylation", cex.lab=0.7)

for (i in seq(along=filteredObjs)) {
    points(start(regionBeta), mcols(regionBeta)[[names(filteredObjs[i])]], col=ifelse(pd$status[i]=="N",1,2), cex=.3, pch=19)
}

for (i in 1:2) {
    curve(predict(fits[[i]], x), from=min(start(regionBeta)), to=max(start(regionBeta)), col=i,add=TRUE,lwd=1.3)
}
legend("bottomright", c("Normal","Tumor"), pch=22,pt.bg=1:2)

indexes <- c(1,4,2,5,3,6)
for (i in indexes) {
    plotPatterns(filteredObjs[[i]], plotRegion,
                 xaxt="n", yaxt="n", ylab="", xlab="",
                 bty="n")
#    mtext(side=3,paste(gsub("CAP_","",names(filteredObjs)[i])), cex=.7)
}
dev.off()


load_all(path)
entrStats <- lapply(filteredObjs, getEntropyStats)
jointGR <- disjoin(unlist(do.call("GRangesList", entrStats)))

makeVec <- function(gr, jointGR, slot) {
    olaps <- findOverlaps(jointGR, gr)
    out <- rep(NA,length(jointGR))
    out[queryHits(olaps)] <- mcols(gr)[[slot]][subjectHits(olaps)]
    out
}
entrMat <- sapply(entrStats, makeVec, jointGR, "entr")
nentrMat <- sapply(entrStats, makeVec, jointGR, "normEntr")
giniMat <- sapply(entrStats, makeVec, jointGR, "gini")
maxEntrMat <- sapply(entrStats, makeVec, jointGR, "maxEntr")

colData=DataFrame(pd)
rownames(colData)=pd$dirname
statsSE <- SummarizedExperiment(rowData=jointGR,
                                colData=colData,
                                assays=SimpleList(entropy=entrMat,
                                    normEntropy=nentrMat,
                                    maxEntropy=maxEntrMat,
                                    gini=giniMat))

pdf(file.path(figdir, "entropy.pdf"),width=9,height=3)
mypar(1,3)

for (i in 1:3) {
    mat = assay(statsSE,"normEntropy")[,c(i+3,i)]
    plot(mat,main="normalized entropy" )
    abline(0,1)
}

for (i in 1:3) {
    mat = assay(statsSE,"normEntropy")[,c(i+3,i)]
    d = mat[,2]-mat[,1]
    hist(d,nc=50,xlab=sprintf("%s-%s", colnames(mat)[2], colnames(mat)[1]),main="normalized entropy")
}

for (i in 1:3) {
    mat = assay(statsSE,"maxEntropy")[,c(i+3,i)]
    plot(mat, main="log num patterns")
    abline(0,1)
}

for (i in 1:3) {
    mat = assay(statsSE,"maxEntropy")[,c(i+3,i)]
    d = mat[,2]-mat[,1]
    hist(d,nc=50,xlab=sprintf("%s-%s", colnames(mat)[2], colnames(mat)[1]),main="log num patterns")
}
dev.off()

jointGR <- disjoin(unlist(do.call(GRangesList, betaGR)))
betaMat <- sapply(betaGR, makeVec, jointGR, "estimatedBeta")
normEntrMat <- sapply(entrStats, makeVec, jointGR, "normEntr")
entrMat <- sapply(entrStats, makeVec, jointGR, "entr")
cpgSE <- SummarizedExperiment(rowData=jointGR,
                              colData=colData,
                              assays=SimpleList(beta=betaMat,
                                  entropy=entrMat,
                                  normEntropy=normEntrMat))

pdf(file.path(figdir,"methEntropy.pdf"),width=9,height=3)
mypar(1,3)

for (i in 1:3) {
    mat = assay(cpgSE, "entropy")[,c(i+3,i)]
    plot(mat,main="entropy")
}

for (i in 1:3) {
    mat = assay(cpgSE,"entropy")[,c(i+3,i)]
    d = mat[,2]-mat[,1]
    hist(d,nc=50,xlab=sprintf("%s-%s", colnames(mat)[2], colnames(mat)[1]),main="entropy")
}

for (i in 1:3) {
    mat = assay(cpgSE, "beta")[,c(i+3,i)]
    plot(mat,main="methylation")
}

for (i in 1:3) {
    mat = assay(cpgSE,"beta")[,c(i+3,i)]
    d = mat[,2]-mat[,1]
    hist(d,nc=50,xlab=sprintf("%s-%s", colnames(mat)[2], colnames(mat)[1]),main="methylation")
}

for (i in 1:6) {
    plot(assay(cpgSE,"beta")[,i], assay(cpgSE,"entropy")[,i],main=pd$dirname[i])
}
dev.off()






################

compflow <-  tapply(patterns(obj)$abundance, patterns(obj)$cid, sum)
patternp <- patterns(obj)$abundance / compflow[as.character(patterns(obj)$cid)]

tmp <- lapply(4:6, function(subj) {
    nindex <- match(paste("CAP","N",subj,sep="_"), names(objs))
    tindex <- match(paste("CAP","T",subj,sep="_"), names(objs))
    nobj <- objs[[nindex]]
    tobj <- objs[[tindex]]
    
    ncomps <- components(nobj)
    tcomps <- components(tobj)
    
    d <- distanceToNearest(ncomps, tcomps)
    keep <- mcols(d)$distance == 0 & width(ncomps)>100 & width(tcomps[subjectHits(d),]) > 100
    iind <- which(keep)
    jind <- subjectHits(d)[keep]
    
    list(comps=ncomps[keep,],
         npats=cbind(npats[[nindex]][iind], npats[[tindex]][jind]),
         ent=cbind(compent[[nindex]][iind], compent[[tindex]][jind]),
         meth=cbind(compMeth[[nindex]][iind], compMeth[[tindex]][jind]))
})

pdf(file.path(figdir,"compcompare.pdf"),height=4,width=12)
mypar(1,3)

  for (i in 1:3) {
    rng <- range(tmp[[i]]$npats, na.rm=TRUE)
    plot(tmp[[i]]$npats, xlim=rng, ylim=rng, xlab="normal num. patterns", ylab="tumor num. patterns",cex.lab=1.7,cex=1.3,pch=19)
    abline(0,1)
    
    rng <- range(tmp[[i]]$ent, na.rm=TRUE)
    plot(tmp[[i]]$ent, xlim=rng, ylim=rng, xlab="normal entropy", ylab="tumor entropy",cex.lab=1.7,cex=1.3,pch=19)
    abline(0,1)
    
    nent <- tmp[[i]]$ent / log(tmp[[i]]$npats)
    rng <- range(nent,na.rm=TRUE)
    plot(nent, xlim=rng, ylim=rng, xlab="normal entropy (normalized)", ylab="tumor entropy (normalized",cex.lab=1.7,cex=1.3,pch=19)
    abline(0,1)
}
dev.off()
  
pdf(file.path(figdir, "compcompare2.pdf"), height=4, width=12)
mypar(1,3)

  for (i in 1:3) {
    nent <- tmp[[i]]$ent / log(tmp[[i]]$npats)
    deltam <- tmp[[i]]$meth[,2]-tmp[[i]]$meth[,1]
    deltae <- nent[,2] - nent[,1]
    plot(deltam, deltae,xlab="methylation difference", ylab="entropy difference",cex.lab=1.7,cex=1.3,pch=19,main=paste("Subject",i))
  }
  dev.off()
  
})




plotRegion2 <- GRanges("chr1",IRanges(start=102600000,end=102700000))
plotComps <- lapply(objs, function(obj) which(components(obj) %over% plotRegion2))

objs2 <- lapply(objs, processMethylpats)

compMeth <- lapply(objs2, componentAvgMeth)

require(epivizr)
mgr <- startEpiviz()

devs <- lapply(seq(along=objs), function(i) mgr$addDevice(components(objs[[i]]), names(objs)[i]))
