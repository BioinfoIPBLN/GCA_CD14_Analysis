## Libraries ####

library("edgeR")
library("genefilter")


## Reading data ####

#reading counts
rawdata<-read.delim(file="ReadCount.tab",header=TRUE,row.names=1,check.names = F)
#reading pheno
targets<-read.table(file = "targets.txt",header = T,stringsAsFactors=F)
#reading gene length
gene.length<-read.table(file="Size.tab",header=T)
idx<-match(rownames(rawdata),gene.length$Gene)
results_counts<-gene.length[idx,]
results_counts[is.na(results_counts$Length),"Length"]<-0

# Assign group

group=targets$Type

## Creating DGE object ####

dge<-DGEList(counts=rawdata, group=group, genes=results_counts) 


## Filtering low expressed genes ####

keep<-rowSums(cpm(data)>1) >= min_group
data<-data[keep,]
data$samples$lib.size <-colSums(data$counts)

## Normalizing ####

dgenorm <- calcNormFactors(dge)
dgenorm <- estimateCommonDisp(dgenorm, robust=TRUE)
dgenorm <- estimateTagwiseDisp(dgenorm)


## Comparison ####

et <- exactTest(dgenorm,pair = comparison)

# Extracting the statistical data order by p-value
top1<-topTags(et, n=nrow(et), adjust.method="BH",sort.by="PValue") 
summary(decideTests(et))
nrow(top1$table[top1$table$FDR<=myFDR & abs(top1$table$logFC)>= myFC,])

write.table(top1$table, file=resultsfile, sep = "\t") 



