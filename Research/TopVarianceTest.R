library(minfi)

workingDir <- "/Users/serraj10/Documents/RCC_Apr26/"
sampleSheet <- file.path(workingDir,"samplesheet.csv")

cmd <- paste("cd", workingDir)
system(cmd)
setwd(workingDir)

targets <- minfi::read.metharray.sheet(base=getwd(), pattern = "samplesheet.csv")
targets$Basename <- file.path(getwd(),targets$SentrixID_Pos)

RGSet <- minfi::read.metharray.exp(base=getwd(), targets=targets, verbose = T, force = T)

detP <- detectionP(RGSet)
colnames(detP) <- RGSet@colData@listData[["Sample_Name"]]
keep <- colMeans(detP) < 0.05
RGSet <- RGSet[, keep]
targets <- targets[keep,]
detP <- detP[, keep]
dropping <- table(keep)["FALSE"]>0
if(!is.na(dropping)){if(dropping==TRUE){message("Dropping probes: ");table(keep)}}
mSetSq <- suppressWarnings(minfi::preprocessQuantile(RGSet))
detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
gset.funnorm <- addSnpInfo(mSetSq[keep, ])
gset.funnorm <- dropLociWithSnps(gset.funnorm,snps=c("SBE","CpG"),maf = 0)#drop loci with snps
annot = getAnnotation(gset.funnorm) #getting annotation file
sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")] # dropping sex probes
gset.funnorm = gset.funnorm[!(rownames(gset.funnorm) %in% sex_probes),]
betas <- minfi::getBeta(gset.funnorm)
colnames(betas) <- RGSet@colData@listData[["Sample_Name"]]

takeTopVariance <- function(betas, topVar){
    var_probes <- apply(betas, 1, var)
    select_var <- names(sort(var_probes, decreasing = T))[topVar]
    top_var_beta <- betas[select_var, ]
    return(top_var_beta)
}

top_var_beta <- takeTopVariance(betas, topVar=1:10000)

annot <- minfi::getAnnotation(RGSet)
topProbes <- rownames(top_var_beta[1:1000,])

probesAnnot <- annot[topProbes,]
print(unique(probesAnnot$chr))

write.csv(probesAnnot, "top_1000_probes_annot.csv")
