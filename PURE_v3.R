setwd("~/PURE-main")
options("scipen"=100, "digits"=4)

# Installing packages ------------------------------------------------------
# If the packages have already been installed, please skip this section
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")
install.packages("ggplot2")
install.packages("plotrix")
install.packages("tidyverse")
install.packages("readxl")
install.packages("xlsx")


# Loading packages --------------------------------------------------------

library(fgsea)
library(ggplot2)
library("plotrix")
library(tidyverse)
library("readxl")
# library("xlsx")

# Benchmarking datasets ---------------------------------------------------


datasets = c("exp1_GSE26487", "exp2_GSE49804", "exp3_GSE86837", "exp4_GSE58434",
                          "exp5_GSE58434", "exp6_GSE11352", "exp7_GSE11352", "exp8_GSE11352", 
                          "exp9_GSE74000", "exp10_GSE12446", "exp11_GSE67266", "exp12_GSE67266", 
                          "exp13_GSE51213", "exp14_GSE58875", "exp15_GSE147507", "exp16_GSE147507")

targetChemical = c("Dexamethasone", "Dexamethasone", "Diethylhexyl Phthalate", "Calcitriol",
                   "Calcitriol", "Estradiol", "Estradiol", "Estradiol", "Acetaminophen",
                   "Estradiol", "Etoposide", "Etoposide", "Dexamethasone", "Copper deficiency", "Methylprednisolone", "Methylprednisolone")


GEOID = c("GSE26487", "GSE49804", "GSE86837", "GSE58434", "GSE58434", "GSE11352", "GSE11352", "GSE11352", 
          "GSE74000", "GSE12446", "GSE67266", "GSE67266", "GSE51213", "GSE58875", "GSE147507", "GSE147507")

organism = c("Human", "Mouse", "Mouse", "Human", "Human", "Human", "Human", "Human", "Human", "Human", "Mouse",  
             "Mouse", "Mouse", "Rat", "Human", "Human")

Hypothesis = c(rep("H1", 14), rep("H2", 2))

# Defining helping functions to load datasets and knowledge base -------------------

create_geneset <- function(network)
{
  tmp <- network %>% group_split(ChemicalID)
  c_name <- sapply(tmp,function(x) x$ChemicalID[1])
  geneset <- lapply(tmp, function(x) unique(x$GeneID))
  names(geneset) <- c_name
  geneset
}

getInputData <- function(identifier, cacheFolder="./data", suffix){
  analysisDataFile <- paste(identifier, suffix, sep='_')
  file = file.path(cacheFolder, analysisDataFile)
  inputData <- read.csv(file=file, header=TRUE,  stringsAsFactors = FALSE)
  return(inputData)
}


# Running PURE algorithm ------------------------------------------------------

run_PURE <- function(gene_sets, gene_set_names, network = network, seed=1, Hypothesis = "H1") {
  set.seed(seed)
  
  res <- lapply(gene_set_names, function(gsn){
    gs <- gene_sets[[gsn]]
    
    DEGenes <- unique(network$GeneID[network$sign == 1 & network$ChemicalID == gsn])
    
    
    l = length(DEGenes)
    k = length(gs) - l
    m = length(unique(network$GeneID[network$sign == 1])) - l
    n = length(unique(network$GeneID[network$sign != 1])) - k
    
    if(Hypothesis == "H1") {  
      p = fisher.test(matrix(c(l,k,m,n), nrow = 2), alternative = "greater")$p.value
    } else {
      p = fisher.test(matrix(c(k,l,n,m), nrow = 2), alternative = "greater")$p.value
    }
    p
  }) %>% unlist() %>% data.frame(stringsAsFactors = F)
  
  colnames(res) <- 'p.value'
  res$padj <- p.adjust(res$p.value, "fdr")
  res$CDT <- gene_set_names
  res <- res %>% tidyr::drop_na()
  res
}

for (i in 1:length(datasets)) {
  set.seed(1)
  chem <- targetChemical[i]
  org <- organism[i]
  if(org == "Rat") org <- "Rattus"
  
  network <- readRDS(paste0("CTD_",org,"_exp.rds"))
  network <- network[,c("ChemicalID", "GeneID", "edge_type")]
  
  all_gene <- unique(network$GeneID)
  
  data <- getInputData(identifier =  datasets[i], suffix = "ga_de.csv")
  data$entrez <- as.character(data$entrez)
  data <- data[data$entrez %in% all_gene, ]
  
  gene_fc <- data$logfc
  names(gene_fc) <- data$entrez
  
  network <- network[network$GeneID %in% data$entrez, ]
  network$sign <- sign(gene_fc[network$GeneID])*network$edge_type

  geneset <- create_geneset(network)

  res <- run_PURE(geneset, names(geneset), network, Hypothesis = Hypothesis[i])
  
  tmp <- readRDS(paste0("CTD_", org, "_exp.rds"))
  tmp <- tmp[, 1:2]
  tmp <- tmp[!duplicated(tmp),]
  rownames(tmp) <- tmp$ChemicalID
  
  res$CDT <- tmp[res$CDT,]$ChemicalName
  
  p_value <- res$padj
  names(p_value) <- res$CDT
  print(rank(p_value, ties.method =  "min")[chem])
  print(sum(p_value < 0.01))
  
  saveRDS(res, file = paste0("./PURE_", datasets[i], "_vOrg.rds"))
}


# Running other benchmarking methods --------------------------------------

run_KS_W <- function(method = "ks", geneSet, FC, gene_p, seed = 1){
  set.seed(seed)
  
  test <- if (method == "ks") ks.test else wilcox.test
  idx <- which(sapply(geneSet, function(gs) sum(gene_p[gs] < 0.05) > 5))
  geneSet <- geneSet[idx]
  
  FC <- FC[unique(unlist(geneSet))]
  
  allGenes <- names(FC)
  
  res <- lapply(geneSet, function(gs){
    
    DEhit <- FC[allGenes[allGenes %in% gs]]
    DEmiss <- FC[allGenes[!allGenes %in% gs]]
    
    if (length(DEhit) == 0 | length(DEmiss) == 0) return(NA)
    
    test(DEhit, DEmiss)$p.value
  }) %>% unlist() %>% data.frame(CDT=names(geneSet), stringsAsFactors = F)
  
  colnames(res) <- c("p.value", "CDT")
  res <- res %>% drop_na()
  res$padj <- p.adjust(res$p.value, "fdr")
  
  res
}

run_ORA <- function(geneSet, gene_p, seed = 1)
{
  set.seed(seed)
  DEGenes <- names(gene_p[gene_p < 0.05])
  
  res <- lapply(geneSet, function(gs){
    
    wBallDraw <- intersect(gs, DEGenes) %>%  length() - 1
    if (wBallDraw < 0) return(1)
    
    wBall <- length(DEGenes)
    bBall <- length(gene_p) - length(DEGenes)
    ballDraw <- length(intersect(gs, names(gene_p)))
    
    1 - phyper(wBallDraw, wBall, bBall, ballDraw)
    
  }) %>% unlist() %>% data.frame(stringsAsFactors = F)
  
  colnames(res) <- 'p.value'
  res$CDT = rownames(res)
  res <- res %>% tidyr::drop_na()
  res$padj <- p.adjust(res$p.value, "fdr")
  
  res
}

for (method in c("ORA", "KS", "Wilcox", "FGSEA")) { 
  for (i in 1:length(datasets)) {
    set.seed(1)
    chem <- targetChemical[i]
    org <- organism[i]
    if(org == "Rat") org <- "Rattus"
    
    network <- readRDS(paste0("CTD_",org,"_exp.rds"))
    # if(org == "Mouse") 
    network <- network[,c("ChemicalID", "GeneID", "edge_type")]
    
    all_gene <- unique(network$GeneID)
    
    if (method == "FGSEA") {
      data <- getInputData(identifier =  datasets[i], suffix = "ga_de.csv")
    } else {
      data <- getInputData(identifier =  datasets[i], suffix = "ga_all.csv")
    }
    
    data$entrez <- as.character(data$entrez)
    data <- data[data$entrez %in% all_gene, ]
    
    gene_fc <- data$logfc
    gene_p <- data$adjpv
    names(gene_fc) <- names(gene_p) <- data$entrez
    
    network <- network[network$GeneID %in% data$entrez, ]
    
    geneset <- create_geneset(network)
    
    if (method == "FGSEA") {
      res <- fgsea(pathways = geneset, stats = gene_fc, nperm = 1e4, minSize = 15)
      res <- data.frame(CDT = res$pathway, p.value = res$pval, padj = res$padj, stringsAsFactors = F)
    }
    if (method == "KS") {
      res <- run_KS_W("ks", geneset, gene_fc, gene_p)
    }
    if (method == "Wilcox") {
      res <- run_KS_W("wilcox", geneset, gene_fc, gene_p)
    }
    if (method == "ORA") {
      res <- run_ORA(geneset, gene_p)
    }
    
    tmp <- readRDS(paste0("CTD_", org, "_exp.rds"))
    tmp <- tmp[, 1:2]
    tmp <- tmp[!duplicated(tmp),]
    rownames(tmp) <- tmp$ChemicalID
    
    res$CDT <- tmp[res$CDT,]$ChemicalName
    
    saveRDS(res, file = paste0("./", method, "_", datasets[i], "_vOrg.rds"))
  }
}


# Collecting the results ---------------------------------------------

methods <- c("PURE", "ORA", "KS", "Wilcox", "FGSEA", "IPA", "IPA_CDT")
result_table <- matrix(nrow = length(datasets), ncol = length(methods)*3) 
colnames(result_table) <- c(paste0("Rank_", methods), paste0("NrSignDrug_", methods), paste0("pAdj_", methods))
rownames(result_table) <- datasets

for (method in methods) {
  for (i in 1:length(datasets)) {
    print(paste0(method, "_", i))
    set.seed(1)
    chem <- targetChemical[i]
    if(!grepl("IPA", method)) {
      res <- readRDS(paste0("./", method, "_", datasets[i], "_vOrg.rds"))
      p_value <- res$padj
      names(p_value) <- res$CDT
      result_table[datasets[i], paste0("Rank_",method)] <- rank(p_value, ties.method = c("average"))[chem]
      result_table[datasets[i], paste0("NrSignDrug_",method)] <- sum(p_value< 0.05)
      result_table[datasets[i], paste0("pAdj_",method)] <- res[res$CDT == chem, "padj"]
    } else {
      res <-  as.data.frame(read_xls(paste0("IPA_analysis/IPA_", datasets[i], ".xls"), sheet = "Sheet1", guess_max = 21474836))
      res <- res[!is.na(res$`Activation z-score`),] 
      if(method == "IPA_CDT") {
        res <- res[grepl("chemical|drug", res$`Molecule Type`),] # for IPA_CDT, filter only chemicals/drugs in molecule type
      }
      res <- res[,c("Upstream Regulator", "Activation z-score", "p-value of overlap")]
      colnames(res) <- c("CDT", "z_score", "p_overlap")
      res$CDT <- toupper(res$CDT)
      z_score <- res$z_score
      names(z_score) <- res$CDT
      
      chem <- case_when(i == 3 ~ "di(2-ethylhexyl) phthalate",
                        i %in% c(6, 7, 8, 10) ~ "beta-estradiol",
                        TRUE ~ chem)
      
      if (toupper(chem) %in% names(z_score)) {
        if(i == 15 | i ==16){
          result_table[datasets[i], paste0("Rank_",method)] <- rank(z_score, ties.method = c("average"))[toupper(chem)]
        } else {
          result_table[datasets[i], paste0("Rank_",method)] <- rank(-z_score, ties.method = c("average"))[toupper(chem)]
        }
        result_table[datasets[i], paste0("pAdj_",method)] <- res[res$CDT == toupper(chem), "p_overlap"]
      }
      if(i == 15 | i ==16){
        result_table[datasets[i], paste0("NrSignDrug_",method)] <- sum(z_score <= -2)
      } else {
        result_table[datasets[i], paste0("NrSignDrug_",method)] <- sum(z_score >= 2)
      }
    }
  }
}

write.csv(result_table, file = "Result_Table.csv", row.names = T)
#Summary for H1
avg <- round(apply(result_table[c(1:14),], 2, mean, na.rm = T),1)
med <- round(apply(result_table[c(1:14),], 2, median, na.rm = T),1)
sd <-round(apply(result_table[c(1:14),], 2, sd, na.rm = T),1)
summary <- paste0(avg,"+/-" ,sd)
names(summary) <- c(paste0("Rank_",methods),paste0("Nr_Sig_",methods), paste0("pv",methods))
print(summary)
#Summary for H2
round(apply(result_table[c(15,16),], 2, mean, na.rm = T),1)

PURE_r = result_table[,"Rank_PURE"]

index = is.na(IPA_CDT_r)
IPA_CDT_r[index] = result_table[index, "NrSignDrug_IPA_CDT"] + 1
IPA_r[index] = result_table[index, "NrSignDrug_IPA"] + 1

wilcox.test(PURE_r, IPA_CDT_r, alternative = "less")
wilcox.test(PURE_r, IPA_r, alternative = "less")
wilcox.test(result_table[,"Rank_PURE"], result_table[,"Rank_ORA"], alternative = "less")
wilcox.test(result_table[,"Rank_PURE"], result_table[,"Rank_KS"], alternative = "less")
wilcox.test(result_table[,"Rank_PURE"], result_table[,"Rank_Wilcox"], alternative = "less")
wilcox.test(result_table[,"Rank_PURE"], result_table[,"Rank_FGSEA"], alternative = "less")

wilcox.test(result_table[,"NrSignDrug_PURE"], result_table[,"NrSignDrug_IPA"], alternative = "less")
wilcox.test(result_table[,"NrSignDrug_PURE"], result_table[,"NrSignDrug_IPA_CDT"], alternative = "less")
wilcox.test(result_table[,"NrSignDrug_PURE"], result_table[,"NrSignDrug_ORA"], alternative = "less")
wilcox.test(result_table[,"NrSignDrug_PURE"], result_table[,"NrSignDrug_KS"], alternative = "less")
wilcox.test(result_table[,"NrSignDrug_PURE"], result_table[,"NrSignDrug_FGSEA"], alternative = "less")
wilcox.test(result_table[,"NrSignDrug_PURE"], result_table[,"NrSignDrug_Wilcox"], alternative = "less")

rankResult <- list(PURE = result_table[,"Rank_PURE"], 
                   ORA = result_table[,"Rank_ORA"], 
                   KS = result_table[,"Rank_KS"], 
                   Wilcoxon = result_table[,"Rank_Wilcox"], 
                   FGSEA = result_table[,"Rank_FGSEA"], 
                   IPA = result_table[,"Rank_IPA"], 
                   `IPA-CDT` = result_table[,"Rank_IPA_CDT"])
rankResult <- lapply(rankResult, FUN = log) 

pdf("RankTarget_vOrg_log_v3.pdf", width = 8, height = 8)
gap.boxplot(rankResult, gap=list(top=c(NA,NA),bottom=c(NA,NA)),
            col=c("firebrick3", "darkorchid1", "cornflowerblue", "chartreuse3", "orange", "#69b3a2", "#488c7c"))
axis(2, labels=c(0), at=c(0))
title(ylab="log(ranks) of the true CDTs", xlab= "Method", col="black")
dev.off()

# Nr of signficant drugs reported

fpResult <- list(PURE = result_table[,"NrSignDrug_PURE"], 
                 ORA = result_table[,"NrSignDrug_ORA"], 
                 KS = result_table[,"NrSignDrug_KS"], 
                 Wilcoxon = result_table[,"NrSignDrug_Wilcox"], 
                 FGSEA = result_table[,"NrSignDrug_FGSEA"], 
                 IPA = result_table[,"NrSignDrug_IPA"], 
                 `IPA-CDT` = result_table[,"NrSignDrug_IPA_CDT"])
pdf("NrSigChem_vOrg_v3.pdf", width = 8, height = 8)

gap.boxplot(fpResult, gap=list(top=c(NA,NA),bottom=c(NA,NA)),
            col=c("firebrick3", "darkorchid1", "cornflowerblue", "chartreuse3", "orange", "#69b3a2", "#488c7c"))
axis(2, labels=c(0), at=c(0))
title(ylab="Number of significant CDTs reported", xlab= "Method", col="black")

dev.off()

