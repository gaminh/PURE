library(fgsea)
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
  analysisDataFile <- paste("upstream", identifier, suffix, sep='_')
  file = file.path(cacheFolder, analysisDataFile)
  inputData <- read.csv(file=file, header=TRUE,  stringsAsFactors = FALSE)
  return(inputData)
}

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
  }) %>% unlist() %>% data.frame(pathway=names(geneSet), stringsAsFactors = F)
  
  colnames(res) <- c("p.value", "pathway")
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
  res$pathway = rownames(res)
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
      res <- data.frame(pathway = res$pathway, p.value = res$pval, padj = res$padj, stringsAsFactors = F)
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
    
    res$pathway <- tmp[res$pathway,]$ChemicalName
    
    saveRDS(res, file = paste0("./result/", method, "_", datasets[i], "_vOrg.rds"))
  }
}