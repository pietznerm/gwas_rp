#!/usr/bin/env Rscript

## Author Dr. Summaira Yasmeen
## Requires LCV (Latent Causal Variable model)  https://github.com/lukejoconnor/LCV/blob/master/README.md

## Scaling LCV analysis on a set of phenoytpes
require(data.table)
require(R.utils)
require(dplyr)

## write a custom function to use LCV method
LCV.func <- function(nameT1, nameT2, ld.file.df=ld.file.df )
{
  ## nameT1 is GWAS for trait 1
  ## nameT2 is GWAS for trait 2
  #Load trait 1 data and calculate Zs
  d1 = na.omit(read.table(gzfile(nameT1),header=TRUE,sep="\t",stringsAsFactors = FALSE))
  
  #Load trait 2 data and calculate Zs
  d2 = na.omit(read.table(gzfile(nameT2),header=TRUE,sep="\t",stringsAsFactors = FALSE))
    
  #Merge
  m = merge(ld.file.df,d1,by="SNP")
  data = merge(m,d2,by="SNP")
  
  #Sort by position 
  data = data[order(data$CHR, data$BP),]
  
  #Flip sign of one z-score if opposite alleles-shouldn't occur with UKB data
  #If not using munged data, will have to check that alleles match-not just whether they're opposite A1/A2
  mismatch = which(data$A1.x!=data$A1.y,arr.ind=TRUE)
  data[mismatch,]$Z.y = data[mismatch,]$Z.y*-1
  data[mismatch,]$A1.y = data[mismatch,]$A1.x
  data[mismatch,]$A2.y = data[mismatch,]$A2.x
  

  #Run LCV
  
  source("..../LCV/R/RunLCV.R")
  
  x = RunLCV(data$L2,data$Z.x,data$Z.y)
  
  T1 <- as.character(strsplit(nameT1, '.sumstats.gz'))
  T2 <- as.character(strsplit(nameT2, '.sumstats.gz'))
  
  new_df <- cbind.data.frame(T1, T2, x$zscore, x$pval.gcpzero.2tailed, x$gcp.pm, x$gcp.pse, x$rho.est, x$rho.err, x$pval.fullycausal[1], x$pval.fullycausal[2], x$h2.zscore[1], x$h2.zscore[2])
  colnames(new_df) <- c('T1', 'T2', 'zscore', 'pval.gcpzero.2tailed','gcp.pm','gcp.pse' , 'rho.est', 'rho.err', 'pval.fullycausal_T1', 'pval.fullycausal_T2', 'h2.zscore_T1', 'h2.zscore_T2')
  setwd('/sc-projects/sc-proj-computational-medicine/people/Summaira/01_rg/LCV_analysis/')
  file.name=paste0(T1,'_', T2, '.lcv')
  write.table(new_df, file = file.name , col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}

## combine chromosome  LD Scores for EUR, into one file. 
## We used scores for EUR population from here https://alkesgroup.broadinstitute.org/LDSCORE/eur_w_ld_chr.tar.bz2
path.ld <- c('..../ldsc/eur_w_ld_chr/')
ld.files <- list.files(path = path.ld , pattern = '.gz')
ld.files.df <- lapply(ld.files, fread)
ld.files.df <- Reduce(full_join,ld.files.df)

## call the function
LCV.func(nameT1 = c('raynaud_gwas_sumstats.gz'), nameT2 = c('date_443.1.sumstats.gz'), ld.file.df = ld.file.df)











