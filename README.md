# PURE
## P(redicting) U(pstream) RE(gulators)

Title: "A novel approach for Predicting Upstream REgulators (PURE) affecting gene expression profile" \
Authors: Tuan-Minh Nguyen<sup>1</sup>, Douglas B. Craig<sup>1,2</sup>, Duc Tran<sup>3</sup>, Tin Nguyen<sup>3</sup>, and Sorin Draghici<sup>1,4,5,\*</sup> 
1. Department of Computer Science, Wayne State University, Detroit, 48202, USA 
2. Department of Oncology, School of Medicine, Wayne State University, Detroit, MI 48201, USA 
3. Department of Computer Science and Engineering, University of Nevada Reno, Reno, 89557, USA 
4. Advaita Bioinformatics, Ann Arbor, MI 48105, USA 
5. Computer and Information Science and Engineering, National Science Foundation, Alexandria, VA 22314, USA \
\* email: sdraghic@nsf.gov 

## 1. Introduction
External factors such as exposure to a chemical, drug, or toxicant (CDT), or conversely, the lack of certain chemicals can cause many diseases. The ability to identify such causal CDTs based on changes in the gene expression profile is extremely important in many studies. Furthermore, the ability to correctly infer CDTs that can revert the gene expression changes induced by a given disease phenotype is a crucial step in drug repurposing. 

We present an approach for Predicting Upstream REgulators (PURE) designed to tackle this challenge. PURE can correctly infer a CDT from the measured expression changes in a given phenotype, as well as correctly identify drugs that could revert disease-induced gene expression changes.

This script, "PURE_v3.R” implements the PURE algorithm and also compare with four classical approaches as well as with the causal analysis used in Ingenuity Pathway Analysis (IPA) on 16 data sets (1 rat, 5 mouse, and 10 human data sets), involving 8 chemicals or drugs.

## 2. Configuration
This script was created using RStudio version 2023.06.0+421 (2023.06.0+421), R version 4.3.0 (2023-04-21), under macOS Monterey 12.0.1. The following packages and versions are needed for running the script:

- readxl, Version 1.4.2
- lubridate, Version 1.9.2
- forcats, Version 1.0.0
- stringr, Version 1.5.0
- dplyr, Version 1.1.2
- purrr, Version 1.0.1
- readr, Version 2.1.4
- tidyr, Version 1.3.0
- tibble, Version 3.2.1
- tidyverse, Version 2.0.0    
- plotrix, Version 3.8-2
- ggplot2, Version 3.4.2
- fgsea, Version 1.26.0
- BiocManager, Version 1.30.21

## 3. Inputs
- **The knowledge base** We used The Comparative Toxicogenomics Database of three species, namely human, mouse, and rat for this study. The knowledge base contains the information about the interactions between CDTs and downstream genes.  It was downloaded from The Comparative Toxicogenomics Database website (https://ctdbase.org)
- **The gene expression profiles** Each data consists of 4 columns: the log(fold-change) of condition’s gene expression compared to control’s ones (column name: logfc), the FDR-adjusted p values of the changes in gene expressions between the condition and control samples (column name: adjpv), entrez ID of the gene (column name: entrez), and the gene symbol (column name: symbol). While some methods, including PURE, only require the information of the DE genes (files with suffix “_ga_de”), others need the information of all the genes measured (files with suffix “_ga_all”). All the gene expression profiles from 16 experiments included in the study are saved in the “data” folder.
- **IPA’s results** The results from IPA were obtained directly from IPA software with default settings and stored in "IPA_analysis_v3" folder.

## 4. Structure
The script consists of 7 sections:
- Installing the packages. 
- Loading packages. 
- Benchmarking datasets. The vectors that contain the information about the data sets are defined in this section, including the dataset name, target CDT, GEO ID, organism, and hypothesis testing for each experiment.
- Defining helping functions to load datasets and knowledge base. The function “create_geneset” creates a list of downstream genes for each CDT for a given knowledge base. The function “getInputData” loads the gene expression profiles into the working environment.
- PURE function (run_PURE). The function takes the following inputs: the sets of downstream genes of the CDTs in the knowledge base (using the “create_geneset” function), the name of the CDTs, the knowledge base, and the hypothesis testing which is either “H1” (default) or H2. The function will return p values and FDR-adjusted p values of the CDTs for the corresponding hypothesis. 
- Implementing other benchmarking methods. The function “run_KS_W” implements the Kolmogorov-Smirnov test (KS) and Wilcoxon test. The parameter “method” takes either the value “ks” for KS or “wilcox” for Wilcoxon test. The other parameters for “run_KS_W” functions are the list of downstream genes for each CDT, the fold-changes, and the p values of the genes. The function “run_ORA” implements the ORA approach and takes the downstream gene sets, and p values of the genes as the parameters. For FGSEA we used the function “fgsea” from the fgsea package. We did not implemented IPA but obtained its results directly from IPA software.
- Collecting the results. The results from running each algorithm are saved in the format of “<method>_<datasets>_vOrg.rds". All the results will be saved in the variable “result_table” which is then saved as "Result_Table.csv" in the same folder. Some statistics mentions in the manuscript such as the average values, the standard deviations, the p value when comparing the performance of PURE and other methods are also calculated in this section. Lastly, we included the code that generates the figures in the manuscript.

## 5. How to run the script
In the terminal, run the following command at the folder containing the script "PURE_v3.R":
```R
Rscript PURE_v3.R
```
The script will prompt user to enter the path where the git repo was cloned.

---
Contact: Please contact tuanminh22@gmail.com if you have any question.