datasets = as.character(c(532, 533, 534, 537, 538, 541, 548, 543, 545, 547, 549, 550, 554, 552, 553, 555, 556, 
                          557, 558, 560, 559, 561, 562, 563, 565, 566, 567, 568, 564, 569, 570, 800))

targetChemical = c("Estradiol", "Estradiol", "Estradiol", "Estradiol", "Estradiol", "Atorvastatin", "Dexamethasone",
                   "Acetaminophen", "Omeprazole", "Prednisolone", "Dexamethasone", "Dexamethasone", "Dexamethasone", "Prednisolone", 
                   "Prednisolone", "Dexamethasone", "Acetaminophen", "Acetaminophen","Tobacco Smoke Pollution", 
                   "Tobacco Smoke Pollution", "Tobacco Smoke Pollution", "Choline deficiency", 
                   "Diethylhexyl Phthalate", "Medroxyprogesterone Acetate", "Etoposide","Etoposide", "Etoposide", 
                   "Etoposide", "Dexamethasone", "Calcitriol", "Calcitriol", "Copper deficiency")

GEOID = c("GSE11352", "GSE11352", "GSE11352", "GSE52649", "GSE12446", "GSE24187", "GSE26487", "GSE74000", "GSE77239", NA, "GSE49804",
          "GSE51213", "GSE51213", NA, "GSE21048", "GSE72907", "GSE68065", "GSE40336", "GSE50254", "GSE50254", "GSE50254",
          "GSE111294", "GSE86837", "GSE68229", "GSE67266", "GSE67266", "GSE67266", "GSE67266", "GSE52778", "GSE58434",
          "GSE58434", "GSE58875")

organism = c("Human", "Human", "Human", "Mouse", "Human", "Human", "Human", "Human", "Human", "Human", "Mouse", "Mouse", "Mouse", "Mouse",
             "Mouse", "Rat", "Rat", "Rat", "Rat", "Rat", "Rat", "Mouse", "Mouse", "Human", "Mouse", "Mouse", "Mouse", "Mouse", "Human", "Human","Human", "Rat")

bad = c(4, 6, 9, 10, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 27, 29)
omit = bad #c(rattus, bad)
datasets = datasets[-omit]
targetChemical = targetChemical[-omit]
GEOID = GEOID[-omit]
organism = organism[-omit]

# setwd("/data/dtran/DrugDiscovery/")
setwd("/Users/minhnguyen/upstream_analysis/")

library(fgsea)
library(tidyverse)

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
