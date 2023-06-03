# Run PURE_v3.R and Other_v3.R scripts before running this file 

setwd("~/PURE-main")

result_table <- matrix(nrow = length(datasets), ncol = 5*3)
methods <- c("ORA", "KS", "Wilcox", "FGSEA", "PURE")
colnames(result_table) <- c(paste0("Rank_", methods), paste0("NrSignDrug_", methods), paste0("pAdj_", methods))
rownames(result_table) <- datasets

for (method in methods) {
  for (i in 1:length(datasets)) {
    set.seed(1)
    chem <- targetChemical[i]

    if(method == "PURE"){
      res <- readRDS(paste0("./", method, "_", datasets[i], "_vOrg.rds"))
    } else {
      res <- readRDS(paste0("./", method, "_", datasets[i], "_vOrg.rds"))
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
      res <- readRDS(paste0("./", method, "_", datasets[i], "_vOrg.rds"))
    } else {
      res <- readRDS(paste0("./", method, "_", datasets[i], "_vOrg.rds"))
    }


    p_value <- res$padj
    names(p_value) <- res$pathway
    result_table05[datasets[i], paste0("Rank_",method)] <- rank(p_value)[chem]
    result_table05[datasets[i], paste0("NrSignDrug_",method)] <- sum(p_value< 0.05)
    result_table05[datasets[i], paste0("pAdj_",method)] <- res[res$pathway == chem, "padj"]
  }
}

# install.packages("kableExtra")
# devtools::install_github("haozhu233/kableExtra")
library(kableExtra)
library(kable)

kable(result_table, style = "latex")

kbl(result_table, booktabs = T) %>%
  kable_styling(latex_options = c("striped", "scale_down"))

apply(result_table05, 2, mean)

PURE_r <- result_table05[,"Rank_PURE"]
ORA_r <- result_table05[,"Rank_ORA"]
KS_r <- result_table05[,"Rank_KS"]
Wilcoxon_r <- result_table05[,"Rank_Wilcox"]
FGSEA_r <- result_table05[,"Rank_FGSEA"]
IPA_r <- c(1, 1, 1, 3, NA, 1, 1, 13, 153, 7, 12, 14, 5, NA)
IPAall_r <- c(1, 1, 816, 25, 26, 1, 1, 1, NA, 10, 22, 21, 34, NA)

mean(IPA_r, na.rm = T)
mean(IPAall_r, na.rm = T)

rankResult <- list(PURE = PURE_r, ORA = ORA_r, KS = KS_r, Wilcoxon = Wilcoxon_r, FGSEA = FGSEA_r, IPA = IPAall_r, `IPA-CDT` = IPA_r)
rankResult <- lapply(rankResult, FUN = log) 

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
