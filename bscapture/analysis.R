library("devtools")
path ="~/Desktop/project/methylFlow/methylFlowr"
load_all(path)

install_github("rafalib","ririzarr")
library("rafalib")
#load("rafalib")

chr <- lapply(seq(1:2), function(j) {
datadir="~/Dropbox/testing/colon"
pd=data.frame(subject=rep(4:6,2),
  status=c("T", "T", "T","N", "N", "N"))
pd$dirname=sprintf("CAP_%s_%d", pd$status, pd$subject)
objs <- lapply(seq(len=nrow(pd)), function(i) {
  curdir <- file.path(datadir,pd$dirname[i],j)
  tmp <- read.methylflow.dir(curdir, pd$dirname[i],has.header=TRUE)
  print(seqlevels(tmp@regions))
  tmp@regions <- renameSeqlevels(tmp@regions,paste("chr",j,sep=""))
  tmp@patterns <- renameSeqlevels(tmp@patterns, paste("chr",j,sep=""))
  tmp@components <- renameSeqlevels(tmp@components, paste("chr",j,sep=""))
  tmp
})



names(objs) <- pd$dirname

compcoverage <- sapply(objs, function(obj) counts(obj, level="component"))

figdir <- file.path(datadir,"figs",j)

pdf(file.path(figdir,"fragment_length.pdf"),height=4,width=6)
mypar(2,3)

for (i in seq(along=objs)) {
  hist(width(components(objs[[i]]))[compcoverage[[i]]>100], nc=10, main=names(objs)[i],xlab="assembled fragment size")
}
dev.off()

npats <- sapply(objs, function(obj) npatterns(obj))

pdf(file.path(figdir,"patt_by_coverage.pdf"), height=4, width=6)
mypar(2,3)

for (i in seq(along=objs)) {
  keep <- width(components(objs[[i]])) > 100
  plot(compcoverage[[i]][keep], npats[[i]][keep],main=names(objs)[i],xlab="coverage",ylab="num. patterns")
}
dev.off()

objs2 <- lapply(objs, processMethylpats)

compent <- sapply(objs, componentEntropy)
compMeth <- lapply(objs2, componentAvgMeth)

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