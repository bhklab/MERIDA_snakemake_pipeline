library(Biobase)
library(PharmacoGx)
library(ggplot2)
library(GSA) 
library(piano)
library(gsubfn)
library(wCI)
library(SummarizedExperiment)

setwd("/root/capsule/code/")

source("src/getBiModalScore.R")
source("src/getPredictions_LOBICO.R")


# Flags to set whether to run long functions
perform_CL_bimodality <- F # flag for recomputing gene expression bimodality on cell lines.
perform_MCC_bimodalGenesSimilarity <- F # flag for recomputing bimodal genes similarity using MCC metric
perform_CTRPv2_training_BimodalGenes <- F # flag for retraining models based on CTRPv2 bimodal genes - RNAseq




######################################
######################################
# Data for cell lines bimodality

CCLE <- readRDS("/root/capsule/data/CCLE.rds")
CCLE <- CoreGx:::.convertCSetMolecularProfilesToSE(CCLE)
genes_mappings <- PharmacoGx::featureInfo(CCLE,"rnaseq")

############
# Processsing Data

ccleRnaSeq_ALL <- t(assay(summarizeMolecularProfiles(CCLE,mDataType = "rnaseq",fill.missing = F)))
ccleTissues <- CCLE@cell[rownames(ccleRnaSeq_ALL),"tissueid",drop=F]
ibx <- which(ccleTissues$tissueid!="haematopoietic_and_lymphoid_tissue")
ccleTissues_WithoutHAEMATO <- ccleTissues[ibx,,drop=F]
ccleRnaSeq <- ccleRnaSeq_ALL[rownames(ccleTissues_WithoutHAEMATO),]
CCLE_tissues <- data.frame(table(ccleTissues[rownames(ccleRnaSeq),]),stringsAsFactors = F)




if(perform_CL_bimodality){
  
  library(parallel)
  cl <- makeCluster(3)
  
  BiModalScores_ccleRnaSeq_ALL <- parApply(cl = cl,FUN = function(x){source("getBiModalScore.R");getBiModalScore_Updated(x)},MARGIN = 2,X = ccleRnaSeq)
  BiModalScores_ccleRnaSeq <- unlist(lapply(BiModalScores_ccleRnaSeq_ALL, "[[",1))
  BiModalScores_ccleRnaSeq <- BiModalScores_ccleRnaSeq[order(BiModalScores_ccleRnaSeq,decreasing = T)]
  
  stopCluster(cl)
}else{
  BiModalScores_ccleRnaSeq <-readRDS("/data/BiModalScores_ccleRnaSeq.rda")
  BiModalScores_ccleRnaSeq_ALL <-readRDS("/data/BiModalScores_ccleRnaSeq_ALL.rda")
}

protCodingGenes <- genes_mappings[names(BiModalScores_ccleRnaSeq),"gene_type"]
names(protCodingGenes) <- names(BiModalScores_ccleRnaSeq)

protCodingGenes <- protCodingGenes[protCodingGenes=="protein_coding"]
protCodingGenes <- names(protCodingGenes)


percentile <- 0.80
ccle_RS_th <- quantile(BiModalScores_ccleRnaSeq[protCodingGenes],probs = percentile)

dens=density(BiModalScores_ccleRnaSeq[protCodingGenes])

# example of one of the top bimodal genes
ibx <- which(genes_mappings$Symbol == "RAB25")


######################################
# Data for TCGA bimodality


#####################
# Fig 1 - c
#####################

TCGA_bimodal_genes <- readRDS("/data/TCGA_bimodality_scores.rda")

TCGA_threshold <- quantile(TCGA_bimodal_genes[protCodingGenes],probs = 0.8)
ccle_Threshold <- quantile(BiModalScores_ccleRnaSeq[protCodingGenes],probs = 0.8)

ibx <- which(BiModalScores_ccleRnaSeq[names(TCGA_bimodal_genes[protCodingGenes])] >= ccle_Threshold & TCGA_bimodal_genes[protCodingGenes] >= TCGA_threshold  )
df1 <- data.frame(data.frame("TCGA"=TCGA_bimodal_genes[protCodingGenes],"CCLE"=BiModalScores_ccleRnaSeq[protCodingGenes],"group"=(TCGA_bimodal_genes[protCodingGenes]>TCGA_threshold & BiModalScores_ccleRnaSeq[protCodingGenes]>ccle_Threshold)))

finalSetOfGenes <- names(ibx)
names(finalSetOfGenes) <- genes_mappings[finalSetOfGenes,"Symbol"]



################################
################################
################################

# Binarizing the gene expression data based on the bimodaly features.
ccleRnaSeq_final <- ccleRnaSeq[,finalSetOfGenes]
colnames(ccleRnaSeq_final) <- finalSetOfGenes

cutoffs_ccle <- unlist(lapply(colnames(ccleRnaSeq_final), function(x){
  val <- BiModalScores_ccleRnaSeq_ALL[[x]][["mix"]][["m.step"]][["mu"]]
  return(mean(val))
}))
names(cutoffs_ccle) <- colnames(ccleRnaSeq_final)

ccleRnaSeq_final_binary <- getBinaryValues(ccleRnaSeq_final,cutoffs_ccle)

ccleRnaSeq_final_binary <- ccleRnaSeq_final_binary[,colSums(is.na(ccleRnaSeq_final_binary))<nrow(ccleRnaSeq_final_binary)]
#Emily's prep for MERIDA 
# compute AAC values 
aac_df <- PharmacoGx::summarizeSensitivityProfiles(CCLE,"aac_recomputed")
#remove na values 
#replace Erlotinib with the drug that you wish to subset for
aac_df <- data.frame(aac_df["Erlotinib",])
colnames(aac_df) <- c("AAC")
aac_df <- na.omit(aac_df)
aac_df <- subset(aac_df, rownames(aac_df) %in% rownames(ccleRnaSeq_final_binary))
rownames(aac_df) <- gsub(" ", "", rownames(aac_df))
write.table(aac_df,file = "/root/capsule/results/AAC_file.txt",sep = " ")
threshold <- aac_df$AAC > 0.2
threshold <- replace(threshold, threshold==TRUE, 1)
aac_df$threshold <- threshold

resistant_sum = 0
sensitive_sum = 0
for(i in 1:nrow(aac_df)) {
  row <- aac_df[i,]
  if (row$threshold == 1) {
    sensitive_sum =  sensitive_sum + abs(row$AAC - log(0.2))**3
  }
  else{
    resistant_sum =  resistant_sum + abs(row$AAC - log(0.2))**3
  }
  # do stuff with row
}
weights <- c()
for(i in 1:nrow(aac_df)) {
  row <- aac_df[i,]
  if (row$threshold == 1) {
    weights[i] <- abs(row$AAC - log(0.2))**3/(2*sensitive_sum)
  }
  else{
    weights[i] <- abs(row$AAC-log(0.2))**3/(2*resistant_sum)
  }
  # do stuff with row
  
}
copy_binary <- data.frame(ccleRnaSeq_final_binary)
rownames(copy_binary) <- gsub(" ", "", rownames(copy_binary))
copy_binary <- subset(copy_binary, rownames(copy_binary) %in% rownames(aac_df))
copy_binary$w <- weights
copy_binary$r <- threshold
write.table(copy_binary,file = "/root/capsule/results/Input_Matrix.txt",sep = " ")

