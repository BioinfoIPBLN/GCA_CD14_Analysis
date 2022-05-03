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
mdsPlot(M, numPositions = 500000co, sampGroups = pdata$Day_sorter, sampNames = pdata$Sample_Name, legendPos = "topright", main = "Meth_Day_Sorter")

    
    
    
 
    
