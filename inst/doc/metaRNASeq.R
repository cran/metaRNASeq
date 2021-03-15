### R code from vignette source 'metaRNASeq.Rnw'

###################################################
### code chunk number 1: metaRNASeq.Rnw:42-43
###################################################
options(width=60)


###################################################
### code chunk number 2: loadparameters
###################################################
library(metaRNASeq)
data(param)
dim(param)
data(dispFuncs)


###################################################
### code chunk number 3: simulateData
###################################################
set.seed(123)
matsim <- sim.function(param = param, dispFuncs = dispFuncs)
sim.conds <- colnames(matsim)
rownames(matsim) <- paste("tag", 1:dim(matsim)[1],sep="")
dim(matsim)


###################################################
### code chunk number 4: extractindivstudy
###################################################
colnames(matsim)
simstudy1 <- extractfromsim(matsim,"study1")
head(simstudy1$study)
simstudy1$pheno
simstudy2 <- extractfromsim(matsim,"study2")


###################################################
### code chunk number 5: DESeq2.indivanalysis
###################################################
 if (requireNamespace("DESeq2", quietly = TRUE)) {
    dds1 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy1$study,
      colData = simstudy1$pheno,design = ~ condition)
    res1 <- DESeq2::results(DESeq2::DESeq(dds1))
    dds2 <- DESeq2::DESeqDataSetFromMatrix(countData = simstudy2$study, 
      colData = simstudy2$pheno,design = ~ condition)
    res2 <- DESeq2::results(DESeq2::DESeq(dds2))
  }


###################################################
### code chunk number 6: storepvalandFC
###################################################
if (exists("res1") && exists("res2"))
{
  rawpval <- list("pval1"=res1[["pvalue"]],"pval2"=res2[["pvalue"]])
  FC <- list("FC1"=res1[["log2FoldChange"]],"FC2"=res2[["log2FoldChange"]])
} else {
  data(rawpval)
  data(FC)
}


###################################################
### code chunk number 7: storeadjpval
###################################################
if (exists("res1") && exists("res2"))
{
  adjpval <- list("adjpval1"=res1[["padj"]],"adjpval2"=res2[["padj"]])
} else {
  data(adjpval)
}

studies <- c("study1", "study2")
DE <- mapply(adjpval, FUN=function(x) ifelse(x <= 0.05, 1, 0))
colnames(DE)=paste("DE",studies,sep=".")


###################################################
### code chunk number 8: pvalDESeq2hist
###################################################
par(mfrow = c(1,2))
hist(rawpval[[1]], breaks=100, col="grey", main="Study 1", xlab="Raw p-values")
hist(rawpval[[2]], breaks=100, col="grey", main="Study 2", xlab="Raw p-values")


###################################################
### code chunk number 9: filteredPval
###################################################
filtered <- lapply(adjpval, FUN=function(pval) which(is.na(pval)))
rawpval[[1]][filtered[[1]]]=NA
rawpval[[2]][filtered[[2]]]=NA


###################################################
### code chunk number 10: pvalDEhist
###################################################
par(mfrow = c(1,2))
hist(rawpval[[1]], breaks=100, col="grey", main="Study 1", 
  xlab="Raw p-values")
hist(rawpval[[2]], breaks=100, col="grey", main="Study 2", 
  xlab="Raw p-values")


###################################################
### code chunk number 11: pvalfishcomb
###################################################
fishcomb <- fishercomb(rawpval, BHth = 0.05)
hist(fishcomb$rawpval, breaks=100, col="grey", main="Fisher method",
  xlab = "Raw p-values (meta-analysis)")


###################################################
### code chunk number 12: pvalinvnorm
###################################################
invnormcomb <- invnorm(rawpval,nrep=c(8,8), BHth = 0.05)   
hist(invnormcomb$rawpval, breaks=100, col="grey", 
  main="Inverse normal method",
  xlab = "Raw p-values (meta-analysis)")    


###################################################
### code chunk number 13: tabDE
###################################################
DEresults <- data.frame(DE, 
  "DE.fishercomb"=ifelse(fishcomb$adjpval<=0.05,1,0),
  "DE.invnorm"=ifelse(invnormcomb$adjpval<=0.05,1,0))
head(DEresults)


###################################################
### code chunk number 14: checkDESeq2
###################################################
signsFC <- mapply(FC, FUN=function(x) sign(x))
sumsigns <- apply(signsFC,1,sum) 
commonsgnFC <- ifelse(abs(sumsigns)==dim(signsFC)[2], sign(sumsigns),0)  


###################################################
### code chunk number 15: filterconflicts
###################################################
unionDE <- unique(c(fishcomb$DEindices,invnormcomb$DEindices))
FC.selecDE <- data.frame(DEresults[unionDE,],do.call(cbind,FC)[unionDE,],
  signFC=commonsgnFC[unionDE], DE=param$DE[unionDE])
keepDE <- FC.selecDE[which(abs(FC.selecDE$signFC)==1),]
conflictDE <- FC.selecDE[which(FC.selecDE$signFC == 0),]
dim(FC.selecDE)
dim(keepDE)
dim(conflictDE)
head(keepDE)


###################################################
### code chunk number 16: filtercheckcache
###################################################
nbtrueconflicts=as.vector(table(conflictDE$DE)[2])


###################################################
### code chunk number 17: filtercheck
###################################################
table(conflictDE$DE)


###################################################
### code chunk number 18: calcul
###################################################
fishcomb_de <- rownames(keepDE)[which(keepDE[,"DE.fishercomb"]==1)] 
invnorm_de <- rownames(keepDE)[which(keepDE[,"DE.invnorm"]==1)] 
indstudy_de <- list(rownames(keepDE)[which(keepDE[,"DE.study1"]==1)], 
                    rownames(keepDE)[which(keepDE[,"DE.study2"]==1)])

IDD.IRR(fishcomb_de,indstudy_de)
IDD.IRR(invnorm_de ,indstudy_de)


###################################################
### code chunk number 19: calcul2
###################################################
x=IDD.IRR(fishcomb_de,indstudy_de)
y=IDD.IRR(invnorm_de ,indstudy_de)


###################################################
### code chunk number 20: venndiagram
###################################################
 if (require("VennDiagram", quietly = TRUE)) {
  venn.plot<-venn.diagram(x = list(study1=which(keepDE[,"DE.study1"]==1),
                                 study2=which(keepDE[,"DE.study2"]==1),
                                 fisher=which(keepDE[,"DE.fishercomb"]==1),
                                 invnorm=which(keepDE[,"DE.invnorm"]==1)),
                        filename = NULL, col = "black",
                        fill = c("blue", "red", "purple","green"),
                        margin=0.05, alpha = 0.6)
  jpeg("venn_jpeg.jpg");
  grid.draw(venn.plot);
  dev.off();
 } 


###################################################
### code chunk number 21: sessionInfo
###################################################
sessionInfo()


