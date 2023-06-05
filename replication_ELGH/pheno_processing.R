# setup
library(tidyverse)
library(Hmisc)
library(ggplot2)
setwd("/genesandhealth/red/BenJacobs/raynaud_gwas//")


# modify covar file 
covars = read_tsv("/genesandhealth/library-red/genesandhealth/GSAv3EAMD/Jul2021_44k_TOPMED-r2_Imputation_b38/GNH.44190.noEthnicOutliers.covariates.20PCs.withS1QST_diabetes.txt")

covars = covars %>% tidyr::separate(SAMPLEID,sep="_",into = c("IID","FID","other")) %>% 
  dplyr::select(FID,IID,6,9,10,12,13:32)
write_tsv(covars,"./outputs/covars.tsv")

# read in pheno data 
pheno_long = read_csv("./inputs/source_of_report_all_datasets_long_format_25_10_22.csv.gz")

# spread to wide format
phenos = pheno_long %>% 
  pivot_wider(id_cols = PseudoNHSNumber, names_from = ICD10, names_prefix = "ICD10__",values_from = source_of_first_report)

phenos %>%
  dplyr::select(contains("I7")) %>%
  glimpse
colnames(phenos)

# define raynauds using 4-digit ICD10
raynauds = phenos %>% 
  filter(!is.na(ICD10__I730))

# create pheno file 
link_file = read_tsv("/genesandhealth/library-red/genesandhealth/2022_05_12_pseudoNHS_oragene_withmissing_DEIDENTIFIED.txt")
pheno_output = link_file %>% 
  dplyr::select(-3) %>%
  mutate(raynaud = ifelse(PseudoNHS %in% raynauds$PseudoNHSNumber,"1","0")) %>% 
  mutate(IID = as.character(OrageneID)) %>%
  left_join(covars %>% dplyr::select(FID,IID),by="IID") %>%
  filter(!is.na(FID))

pheno_output = pheno_output %>%
  dplyr::select(FID,IID,-PseudoNHS,-OrageneID,3:ncol(pheno_output))
  
write_tsv(pheno_output,"./outputs/phenos.tsv")

table(pheno_output$raynaud)
