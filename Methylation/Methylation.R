## Libraries ####

library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(limma)
library(minfi)
library(minfiData)
library(data.table)

## Reading data ####

baseDir <- "Raw_data"
targets <- read.metharray.sheet(baseDir)

rgset <- minfi::read.metharray.exp(
  targets = targets, verbose = TRUE,
  force = TRUE)

#pheno
pdata <- as.data.frame(pData(rgset))
#SNPs
snps <- data.frame(getSnpBeta(rgset))

#Quality Control before normalization
quality_copntrol <- as.data.frame(minfi::getQC(minfi::preprocessRaw(rgset))) %>%
  dplyr::mutate(
    badSampleCutoff = 10,
    Sample = pdata$Sample_Name,
    Quality = ifelse(
      .data$mMed < badSampleCutoff |
        .data$uMed < badSampleCutoff,
      "Suboptimal",
      "OK"
    )
    
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
    
# remove SNPs and CHs
gset <- minfi::dropLociWithSnps(gset)
gset <- minfi::dropMethylationLoci(gset,
                                   dropRS = TRUE,
                                   dropCH = TRUE)
# remove chromosomes
gset <- gset[rownames(minfi::getAnnotation(gset))[!(minfi::getAnnotation(gset)$chr %in%
                                                      c("chrX", "chrY"))], ]
# remove bad probes
gset <- gset[minfi::featureNames(gset)[!(minfi::featureNames(gset) %in% bad_pos)],]

# Quality control after normalization
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

# Multidimensinal scaling plot - Batch effect

mdsPlot(M, numPositions = 500000, sampGroups = pdata$Group, sampNames = pdata$Sample_Name, legendPos =  "topright", main = "Meth_Disease")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Disease_Group, sampNames = pdata$Sample_Name, legendPos =  "topright", main = "Meth_Disease")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Sex, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Sex")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Age, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Age")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Slide, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Slide")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Array, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Array")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Sample_plate, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Sample_Plate")
mdsPlot(M, numPositions = 500000, sampGroups = pdata$Day_sorter, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Day_Sorter")
   
# PCA
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
    
# Statistical analysis
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


# GCA vs CTRL
GCA_vs_CTRL <- T_limma(input=na.omit(M), pdata = pdata, sample_group =  pdata$Group, age = NULL, arrayweights = T, trend = T, robust = T,  testgroup = "GCA", controlgroup = "CTRL")
GCA_vs_CTRL_FDR <- GCA_vs_CTRL[GCA_vs_CTRL$adj.P.Val < 0.05,] 
GCA_vs_CTRL_Hyper <- GCA_vs_CTRL_FDR[GCA_vs_CTRL_FDR$logFC > 0,] 
GCA_vs_CTRL_Hypo <- GCA_vs_CTRL_FDR[GCA_vs_CTRL_FDR$logFC < 0,]
    
    
    
    
 
    
