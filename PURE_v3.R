options("scipen"=100, "digits"=4)
library(ggplot2)
library("plotrix")
library(tidyverse)

datasets = c("exp1_GSE26487", "exp2_GSE49804", "exp3_GSE86837", "exp4_GSE58434",
                          "exp5_GSE58434", "exp6_GSE11352", "exp7_GSE11352", "exp8_GSE11352", 
                          "exp9_GSE74000", "exp10_GSE12446", "exp11_GSE67266", "exp12_GSE67266", 
                          "exp13_GSE51213", "exp14_GSE58875")

targetChemical = c("Dexamethasone", "Dexamethasone", "Diethylhexyl Phthalate", "Calcitriol",
                   "Calcitriol", "Estradiol", "Estradiol", "Estradiol", "Acetaminophen",
                   "Estradiol", "Etoposide", "Etoposide", "Dexamethasone", "Copper deficiency")


GEOID = c("GSE26487", "GSE49804", "GSE86837", "GSE58434", "GSE58434", "GSE11352", "GSE11352", "GSE11352", 
          "GSE74000", "GSE12446", "GSE67266", "GSE67266", "GSE51213", "GSE58875")

organism = c("Human", "Mouse", "Mouse", "Human", "Human", "Human", "Human", "Human", "Human", "Human", "Mouse",  
             "Mouse", "Mouse", "Rat")

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

run_PURE <- function(gene_sets, gene_set_names, network = network, seed=1) {
  set.seed(seed)
  
  res <- lapply(gene_set_names, function(gsn){
    gs <- gene_sets[[gsn]]
    
    DEGenes <- unique(network$GeneID[network$sign == 1 & network$ChemicalID == gsn])
    
    wBallDraw <- length(DEGenes) - 1
    if (wBallDraw < 0) return(1)
    
    wBall <- length(unique(network$GeneID[network$sign == 1]))
    bBall <- length(unique(network$GeneID[network$sign != 1]))
    ballDraw <- length(gs)
    
    phyper(wBallDraw, wBall, bBall, ballDraw, lower.tail = F)
    
    ## One can use fisher.test as follows:
    # fisher.test(matrix(c(length(DEGenes), ballDraw - length(DEGenes), 
    #                     wBall - length(DEGenes), 
    #                     bBall - ballDraw + length(DEGenes)), nrow = 2), 
    #            alternative = "greater")$p.value
  }) %>% unlist() %>% data.frame(stringsAsFactors = F)
  
  colnames(res) <- 'p.value'
  res$padj <- p.adjust(res$p.value, "fdr")
  res$pathway <- gene_set_names
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

  res <- run_PURE(geneset, names(geneset), network)
  
  tmp <- readRDS(paste0("CTD_", org, "_exp.rds"))
  tmp <- tmp[, 1:2]
  tmp <- tmp[!duplicated(tmp),]
  rownames(tmp) <- tmp$ChemicalID
  
  res$pathway <- tmp[res$pathway,]$ChemicalName
  
  p_value <- res$padj
  names(p_value) <- res$pathway
  print(rank(p_value, ties.method =  "min")[chem])
  print(sum(p_value < 0.01))
  
  saveRDS(res, file = paste0("./PURE_", datasets[i], "_vOrg.rds"))
}