options("scipen"=100, "digits"=4)
library(ggplot2)
library("plotrix")
library(tidyverse)

datasets = as.character(c(548, 549, 562, 569, 570, 532, 533, 534, 543, 538, 568, 566, 554, 800))

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
  
  # load(paste0("PropagatedNetwork_", org, "_v2.RData"))
  network <- readRDS(paste0("CTD_",org,"_exp.rds"))
  # if(org == "Mouse") 
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
  
  # saveRDS(res, file = paste0("./result/PURE_", datasets[i], "_vOrg.rds"))
}


stop("Collecting results. Run other methods first:")

result_table <- matrix(nrow = length(datasets), ncol = 5*3)
methods <- c("ORA", "KS", "Wilcox", "FGSEA", "PURE")
colnames(result_table) <- c(paste0("Rank_", methods), paste0("NrSignDrug_", methods), paste0("pAdj_", methods))
rownames(result_table) <- datasets

for (method in methods) {
  for (i in 1:length(datasets)) {
    set.seed(1)
    chem <- targetChemical[i]

    if(method == "PURE"){
      res <- readRDS(paste0("./result/", method, "_", datasets[i], "_vOrg.rds"))
    } else {
      res <- readRDS(paste0("./result/", method, "_", datasets[i], "_vOrg.rds"))
    }
    
    
    p_value <- res$padj
    names(p_value) <- res$pathway
    result_table[datasets[i], paste0("Rank_",method)] <- rank(p_value)[chem]
    result_table[datasets[i], paste0("NrSignDrug_",method)] <- sum(p_value< 0.01)
    result_table[datasets[i], paste0("pAdj_",method)] <- res[res$pathway == chem, "padj"]
  }
}

result_table05 <- matrix(nrow = length(datasets), ncol = 5*3)
methods <- c("ORA", "KS", "Wilcox", "FGSEA", "PURE")
colnames(result_table05) <- c(paste0("Rank_", methods), paste0("NrSignDrug_", methods), paste0("pAdj_", methods))
rownames(result_table05) <- datasets

for (method in methods) {
  for (i in 1:length(datasets)) {
    set.seed(1)
    chem <- targetChemical[i]
    
    if(method == "PURE"){
      res <- readRDS(paste0("./result/", method, "_", datasets[i], "_vOrg.rds"))
    } else {
      res <- readRDS(paste0("./result/", method, "_", datasets[i], "_vOrg.rds"))
    }
    
    
    p_value <- res$padj
    names(p_value) <- res$pathway
    result_table05[datasets[i], paste0("Rank_",method)] <- rank(p_value)[chem]
    result_table05[datasets[i], paste0("NrSignDrug_",method)] <- sum(p_value< 0.05)
    result_table05[datasets[i], paste0("pAdj_",method)] <- res[res$pathway == chem, "padj"]
  }
}

install.packages("kableExtra")
# devtools::install_github("haozhu233/kableExtra")
library(kableExtra)
library(kable)
kable(result_table, style = "latex")

kbl(result_table, booktabs = T) %>%
  kable_styling(latex_options = c("striped", "scale_down"))

apply(result_table05, 2, mean)
apply(result_table05, 2, median) # \multicolumn{3}{c}{Median}	 &\cellcolor{green}1.5	&	20	&	24	&	53.8	&	27	&	15.5& 4\\ 


PURE_r <- result_table05[,"Rank_PURE"]
ORA_r <- result_table05[,"Rank_ORA"]
KS_r <- result_table05[,"Rank_KS"]
Wilcoxon_r <- result_table05[,"Rank_Wilcox"]
FGSEA_r <- result_table05[,"Rank_FGSEA"]
IPA_r <- c(1, 1, 1, 3, NA, 1, 1, 13, 153, 7, 12, 14, 5, NA)
IPAall_r <- c(1, 1, 816, 25, 26, 1, 1, 1, NA, 10, 22, 21, 34, NA)

median(IPA_r, na.rm = T)
median(IPAall_r, na.rm = T)


PURE_r <- log(PURE_r)
IPA_r <- log(IPA_r)
IPAall_r <- log(IPAall_r)
ORA_r <- log(ORA_r)
KS_r <- log(KS_r)
Wilcoxon_r <- log(Wilcoxon_r)
FGSEA_r <- log(FGSEA_r)

from <- 0
to <- 6

rankResult <- list(PURE = PURE_r, ORA = ORA_r, KS = KS_r, Wilcoxon = Wilcoxon_r, FGSEA = FGSEA_r, IPA = IPAall_r, `IPA-CDT` = IPA_r)

mean(c(1, 1, 816, 25, 26, 1, 1, 1, NA, 10, 22, 21, 34, NA), na.rm = T)

pdf("RankTarget_vOrg_log_v3.pdf", width = 8, height = 8) 
gap.boxplot(rankResult, gap=list(top=c(NA,NA),bottom=c(NA,NA)),
            col=c("firebrick3", "darkorchid1", "cornflowerblue", "chartreuse3", "orange", "#69b3a2", "#488c7c"))
axis(2, labels=c(0), at=c(0))
title(ylab="log(ranks) of the true CDTs", xlab= "Method", col="black")
dev.off()

# Nr of signficant drugs reported
PURE_fp <- result_table05[,"NrSignDrug_PURE"]
IPA_fp <- c(6, 13, 22, 31, 14, 14, 10, 174, 14, 53, 61, 69, 35, 5)
IPAall_fp <- c(23, 18, 191, 169, 198, 79, 125, 324, 30, 226, 155, 341, 531, 16)
ORA_fp <- result_table05[,"NrSignDrug_ORA"]
KS_fp <- result_table05[,"NrSignDrug_KS"]
Wilcoxon_fp <- result_table05[,"NrSignDrug_Wilcox"]
FGSEA_fp <- result_table05[,"NrSignDrug_FGSEA"]

fpResult <- list(PURE = PURE_fp, ORA = ORA_fp, KS = KS_fp, Wilcoxon = Wilcoxon_fp, FGSEA = FGSEA_fp, IPA = IPAall_fp, `IPA-CDT` = IPA_fp)
pdf("NrSigChem_vOrg_v3.pdf", width = 8, height = 8) 

gap.boxplot(fpResult, gap=list(top=c(NA,NA),bottom=c(NA,NA)),
            col=c("firebrick3", "darkorchid1", "cornflowerblue", "chartreuse3", "orange", "#69b3a2", "#488c7c"))
axis(2, labels=c(0), at=c(0))
title(ylab="Number of significant CDTs reported", xlab= "Method", col="black")

dev.off()

lapply(fpResult, mean)
lapply(fpResult, sd)
