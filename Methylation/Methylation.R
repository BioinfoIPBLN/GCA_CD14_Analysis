## Libraries ####

library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
library(limma)
library(minfi)
library(minfiData)
library(data.table)
library(factoextra)

## Reading data ####
baseDir <- "Raw_data"
targets <- read.metharray.sheet(baseDir)
rgset <- minfi::read.metharray.exp(
  targets = targets, verbose = TRUE,
  force = TRUE)

#pheno
pdata <- as.data.frame(pData(rgset))
# SNPs
snps <- data.frame(getSnpBeta(rgset))

#Quality Control before normalization
quality_control <- as.data.frame(minfi::getQC(minfi::preprocessRaw(rgset))) %>%
  dplyr::mutate(
    badSampleCutoff = 10,
    Sample = pdata$Sample_Name,
    Quality = ifelse(
      .data$mMed < badSampleCutoff |
        .data$uMed < badSampleCutoff,
      "Suboptimal",
      "OK"))
qcReport(rgset, sampNames=pdata$Sample_Name, sampGroups=pdata$Disease_Group, pdf="qcReport.pdf")
densityPlot(rgset, sampGroups = pdata$Disease_Group, main = "METHYLATION", legend = TRUE, xlim = c(0,1.0), ylim =c(0,4.5))
 
#Bad probes
detP <- minfi::detectionP(rgset)
bad_pos <- row.names(as.data.frame(detP))[rowMeans(detP) > 0.01]

#Normalization
gset <- minfi::preprocessQuantile(minfi::preprocessNoob(rgset))

#Sex prediction
sex_info <- as.data.frame(minfi::pData(gset)[, c("xMed", "yMed", "predictedSex")])
sex_info$sample <- pdata$Sample_Name
    
#Remove SNPs and CHs
gset <- minfi::dropLociWithSnps(gset)
gset <- minfi::dropMethylationLoci(gset,
                                   dropRS = TRUE,
                                   dropCH = TRUE)
#Remove chromosomes
gset <- gset[rownames(minfi::getAnnotation(gset))[!(minfi::getAnnotation(gset)$chr %in%
                                                      c("chrX", "chrY"))], ]
#Remove bad probes
gset <- gset[minfi::featureNames(gset)[!(minfi::featureNames(gset) %in% bad_pos)],]

#Quality control after normalization
qcReport_gset <- function(input=gset, sampNames = x, path = ""){
  pdf(path)
  densityPlot(input)
  densityBeanPlot(input, sampNames = sampNames)
  dev.off()}
qcReport_gset(input = as.matrix(getBeta(gset)), sampNames = pdata$Sample_Name, path = "QC.pdf")
    
#Get Beta and M values
#Eliminate Inf values from 
Beta <- as.data.frame(getBeta(gset))
Beta[Beta == 0]<-0.0001
Beta[Beta > 0.98999999]<-0.989
Beta <- na.omit(Beta)
M = lumi::beta2m(Beta)

#Multidimensinal scaling plot - Batch effect
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Group, sampNames = pdata$Sample_Name, legendPos =  "topright", main = "Meth_Disease")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Disease_Group, sampNames = pdata$Sample_Name, legendPos =  "topright", main = "Meth_Disease")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Sex, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Sex")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Age, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Age")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Slide, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Slide")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Array, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Array")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Sample_plate, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Sample_Plate")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Day_sorter, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Day_Sorter")
   
#PCA
M1 <- data.frame(t(M))
PCobj = prcomp(M1, retx = T, center = T, scale. = T)
attributes(PCobj) 
PCs = PCobj$x 
summary(PCobj)
autoplot(PCobj, data = pdata, colour = "Group", label.size = 3, shape = FALSE)
autoplot(PCobj, data = pdata, colour = "Disease_Group", label.size = 3, shape = FALSE)
autoplot(PCobj, data = pdata, colour = "Sex", label.size = 3, shape = FALSE)
autoplot(PCobj, data = pdata, colour = "Age", label.size = 3, shape = FALSE)
autoplot(PCobj, data = pdata, colour = "Slide", label.size = 3, shape = FALSE)
autoplot(PCobj, data = pdata, colour = "Array", label.size = 3, shape = FALSE)
autoplot(PCobj, data = pdata, colour = "Sample_plate", label.size = 3, shape = FALSE)
autoplot(PCobj, data = pdata, colour = "Sample_Well", label.size = 3, shape = FALSE)
autoplot(PCobj, data = pdata, colour = "Day_sorter", label.size = 3, shape = FALSE)
    
#Statistical analysis
Function_limma <- function(input = M, pdata = pdata, sample_group = sample_group, sex = NULL, age = NULL, 
                    batch = NULL, arrayweights = F, trend = F, robust = F, testgroup = "GCA", controlgroup = "CTRL"){
  library(limma)
  if(is.null(sex)){
    if(is.null(age)){
      if(is.null(batch)){message("Covariate detected: NULL")
        design = model.matrix(~0+sample_group)
      } else{message("Covariate detected: batch")
        design = model.matrix(~0+sample_group + batch)
      }
    } else{
      message("Covariate detected: age")
      if(is.null(batch)){
        design = model.matrix(~0+sample_group + age)
      } else{message("Covariate detected: batch")
        design = model.matrix(~0+sample_group + batch + age)
      }
    }
  } else{message("Covariate detected: sex")
    if(is.null(age)){
      if(is.null(batch)){
        design = model.matrix(~0+sample_group + sex)
      } else{message("Covariate detected: batch")
        design = model.matrix(~0+sample_group + batch + sex)
      }
    } else{
      if(is.null(batch)){
        design = model.matrix(~0+sample_group + age + sex)
      } else{message("Covariate detected: age")
        design = model.matrix(~0+sample_group + batch + age + sex)
      }
    }
  }
  
  if (!isTRUE(arrayweights)){message("No array weights detected")
    fit <- lmFit(input, design)
    contrast <-paste(paste0("sample_group", testgroup), 
                     paste0("sample_group", controlgroup), sep="-")
    contrast_tissue.matrix <- makeContrasts(contrasts = contrast, levels=design)
    fit_Patients = contrasts.fit(fit, contrast_tissue.matrix)
    fit_Patients <- eBayes(fit_Patients, trend = trend)
    topTable(fit_Patients, coef=1,adjust.method="fdr", n=Inf)
  }else {message("Fitting array weights...")
    array_weight = arrayWeights(input, design=design)
    corweight = asMatrixWeights(array_weight, dim(input))
    fit <- lmFit(input, design,weights=corweight)
    contrast <-paste(paste0("sample_group", testgroup), 
                     paste0("sample_group", controlgroup), sep="-")
    contrast_tissue.matrix <- makeContrasts(contrasts = contrast, levels=design)
    fit_Patients = contrasts.fit(fit, contrast_tissue.matrix)
    fit_Patients <- eBayes(fit_Patients, trend = trend, robust = robust)
    topTable(fit_Patients, coef=1,adjust.method="fdr", n=Inf)
  }
}

#GCA vs CTRL
GCA_vs_CTRL <- Function_limma(input=na.omit(M), pdata = pdata, sample_group =  pdata$Group, age = NULL, arrayweights = T, trend = T, robust = T,  testgroup = "GCA", controlgroup = "CTRL")
GCA_vs_CTRL_FDR <- GCA_vs_CTRL[GCA_vs_CTRL$adj.P.Val < 0.05,] 
GCA_vs_CTRL_Hyper <- GCA_vs_CTRL_FDR[GCA_vs_CTRL_FDR$logFC > 0,] 
GCA_vs_CTRL_Hypo <- GCA_vs_CTRL_FDR[GCA_vs_CTRL_FDR$logFC < 0,]
 
#β-value differentially
Function_merge <- function(x = beta, y = FDR, by = 0){
  matrix <- merge(x, y, by=by)
  rownames(matrix)<-matrix[,1]
  matrix[,-1]
}
Beta_merged <- Function_merge(x = Beta, y =GCA_vs_CTRL , by=0)    
CTRL <- pdata$Sample_Name[pdata$Disease_Group=="CTRL"]
GCA <- pdata$Sample_Name[pdata$Disease_Group=="GCA"]
Beta_merged$GCA_CTRL_dif <- apply(beta_merged[, GCA],1,mean) -apply(beta_merged[,CTRL],1,mean)
Beta_merged$CTRL_mean<-apply(beta_merged[,colnames(beta_merged) %in% CTRL], 1, mean)
Beta_merged$GCA_mean<-apply(beta_merged[,colnames(beta_merged) %in% GCA], 1, mean)
DMPs_BETA <- subset(Beta_merged, (Beta_merged$GCA_CTRL1_adj.P.Val < 0.05))
    
#Annotation 
Annotation <- as.data.frame(getAnnotation(gset))
    
#RefGene_Group
UCSC_RefGene_Group <- data.frame(annotation$UCSC_RefGene_Group, rownames(annotation))
RefGene_Group <- UCSC_RefGene_Group[UCSC_RefGene_Group$rownames.annotation. %in% rownames(GCA_vs_CTRL_FDR),]
Borrar <- c("rownames.annotation.")
RefGene_Group2 <- data.frame(RefGene_Group[, !(names(RefGene_Group)%in% Borrar)]) 
rownames(RefGene_Group2) <- RefGene_Group$rownames.annotation.
GRUPO = data.frame(RefGene_Group2[c(rownames(GCA_vs_CTRL_FDR)),])
rownames(GRUPO) <- rownames(GCA_vs_CTRL_FDR)

#Relation_to_Island
Relation_to_Island <- data.frame(annotation$Relation_to_Island, rownames(annotation))
to_Island <- Relation_to_Island[Relation_to_Island$rownames.annotation. %in% rownames(GCA_vs_CTRL_FDR),] 
Borrar <- c("rownames.annotation.")
to_Island2 <- data.frame(to_Island[, !(names(to_Island)%in% Borrar)])
rownames(to_Island2) <- to_Island$rownames.annotation.
ISLA = data.frame(to_Island2[c(rownames(GCA_vs_CTRL_FDR)),])
rownames(ISLA) <- rownames(GCA_vs_CTRL_FDR)

#Genes
Name_genes <- data.frame(annotation$UCSC_RefGene_Name, rownames(annotation))
Genes <- Name_genes[Name_genes$rownames.annotation. %in% rownames(GCA_vs_CTRL_FDR),]
Borrar <- c("rownames.annotation.")
Genes2 <- data.frame(Genes[, !(names(Genes)%in% Borrar)])
rownames(Genes2) <- Genes$rownames.annotation.
GENES = data.frame(Genes2[c(rownames(GCA_vs_CTRL_FDR)),])
rownames(GENES) <- rownames(GCA_vs_CTRL_FDR)

#Chromosome
Chromosome <- data.frame(annotation$chr, rownames(annotation))
Chromosome2 <- Chromosome[Chromosome$rownames.annotation.. %in% rownames(GCA_vs_CTRL_FDR),]
Borrar <- c("rownames.annotation.")
Chromosome2 <- data.frame(Chromosome[, !(names(Chromosome)%in% borrar)])
rownames(Chromosome2) <- Chromosome$rownames.annotation.
CHR = data.frame(Chromosome2[c(rownames(GCA_vs_CTRL_FDR)),])
rownames(CHR) <- rownames(GCA_vs_CTRL_FDR)

#Position
BP <- data.frame(annotation$pos, rownames(annotation))
BP2 <- BP[BP$rownames.annotation. %in% rownames(GCA_vs_CTRL_FDR),]
Borrar <- c("rownames.annotation.")
BP2 <- data.frame(BP[, !(names(BP)%in% Borrar)])
rownames(BP2) <- BP$rownames.annotation.
BP_F = data.frame(BP2[c(rownames(GCA_vs_CTRL_FDR)),])
rownames(BP_F) <- rownames(GCA_vs_CTRL_FDR)

#Final table
DMPs_BETA_ORDEN_FDR = DMPs_BETA[c(rownames(GCA_vs_CTRL_FDR)),]
GCA_CTRL_FDR_P_FDR <- GCA_vs_CTRL_FDR[,c(1,4,5)]
TABLA_FINAL_GCA_CTRL <- cbind(CHR, BP_F, GENES, GRUPO, ISLA, DMPs_BETA_ORDEN_FDR, GCA_CTRL_FDR_P_FDR)
TABLA_FINAL_GCA_CTRL$METHYLATION = ifelse(TABLA_FINAL_GCA_CTRL$FDR < 0.05 & abs(TABLA_FINAL_GCA_CTRL$Dif_Beta) >= 0, 
                                  ifelse(TABLA_FINAL_GCA_CTRL$Dif_Beta > 0 ,'Hypermethylated','Hypomethylated'),'Stable')
write.table(TABLA_FINAL_GCA_CTRL, file = "TABLA_FINAL_GCA_CTRL.txt", sep = "\t", quote = F)
    

 
    
