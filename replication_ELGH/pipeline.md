###############################
# Step 1 - Pheno & Covars     #
###############################

````unix
Rscript pheno_processing.R
````

###############################
# Step 2 - Gentype QC	      #
###############################
# Prior to this QC extensive QC has been run by Sanger

```unix
cd /genesandhealth/red/BenJacobs/reproductive_gwas/outputs/
rm good_vars_maf_0.01_info_0.7
# filter by INFO > 0.7 and MAF > 0.01
for i in {1..22}
do 
	awk '{if($7>0.01 && $10 > 0.7) print $1,$7,$9,$10}' /genesandhealth/library-red/genesandhealth/GSAv3EAMD/Jul2021_44k_TOPMED-r2_Imputation_b38/topmed-r2_merged_version03/chr$i\_variant_info.tab >> good_vars_maf_0.01_info_0.7
echo "finishied copying list of good vars for chrom $i\ "
done

# convert from BGEN to PLINK1 binary and perform basic QC

# copy big pgen 
cp /genesandhealth/library-red/genesandhealth/GSAv3EAMD/Jul2021_44k_TOPMED-r2_Imputation_b38/topmed-r2_merged_version03/chrALLincX.dose.merged_INFO0.3_MAF0.00001_F_MISSING0.2.8bit.sorted* /home/ivm/

cd /genesandhealth/red/BenJacobs/reproductive_gwas/outputs/
plink2 --pfile /home/ivm/chrALLincX.dose.merged_INFO0.3_MAF0.00001_F_MISSING0.2.8bit.sorted \
--extract good_vars_maf_0.01_info_0.7 \
--out /home/ivm/chr_all_filtered \
--maf 0.01 \
--hwe 1e-10 \
--geno 0.1 \
--mind 0.1 \
--chr 1-22 \
--make-pgen


# convert IIDs
````R
# modify psam file
geno = read_tsv("/home/ivm/chr_all_filtered.psam")
geno = geno %>% 
mutate(old_iid = IID) %>%
mutate(old_fid = 1) %>%
separate(IID,sep="_",into=c("IID","FID","other")) %>% 
select(old_fid,old_iid,FID,IID)
write_tsv(geno,"./outputs/new_iids")
````

# rename IIDs
# rename SNPs as chr:pos
cd /genesandhealth/red/BenJacobs/reproductive_gwas/outputs/
plink2 --pfile /home/ivm/chr_all_filtered \
--set-all-var-ids @:# \
--update-ids new_iids \
--out /home/ivm/post_qc_all_chrs_for_regenie \
--make-pgen



###############################
# Step 2 - GWAS               #
###############################

````unix
cd /genesandhealth/red/BenJacobs/raynaud_gwas/

# copy INFO stats
echo "copying INFO stats"
rm ./outputs/info_1_snps
for i in {1..22}
do 
	awk '{if($7>0.01 && $10 >0.99 1) print $1,$7,$9,$10}' /genesandhealth/library-red/genesandhealth/GSAv3EAMD/Jul2021_44k_TOPMED-r2_Imputation_b38/topmed-r2_merged_version03/chr$i\_variant_info.tab >> ./outputs/info_1_snps
echo "finishied copying list of good vars for chrom $i\ "
done
````

Rename in R
````R
library(tidyverse)
snps = read_table("./outputs/info_1_snps",col_names=F) %>%
  tidyr::separate(X1,sep = ":",into=c("chr","pos","A1","A2"))

snps$CHR = str_remove(snps$chr,"chr")
snps$ID = paste0(snps$CHR,":",snps$pos)
snps = snps %>% dplyr::select(ID)
write_tsv(snps,"./outputs/snps_for_step1_regenie",col_names = F)

````

### REGENIE 
````unix
cd /genesandhealth/red/BenJacobs/raynaud_gwas/

regenie \
--step 1 \
--pgen /home/ivm/post_qc_all_chrs_for_regenie \
--extract ./outputs/snps_for_step1_regenie \
--phenoFile ./outputs/phenos.tsv \
--covarFile ./outputs/covars.tsv \
--bt \
--bsize 1000 \
--lowmem \
--out ./outputs/raynaud_regenie_step1


### step 2
cd /genesandhealth/red/BenJacobs/raynaud_gwas/
regenie \
--step 2 \
--pgen /home/ivm/post_qc_all_chrs_for_regenie \
--phenoFile ./outputs/phenos.tsv \
--covarFile ./outputs/covars.tsv \
--bsize 1000 \
--bt \
--firth --approx \
--pThresh 0.01 \
--pred ./outputs/raynaud_regenie_step1_pred.list \
--out ./outputs/genes_and_health_raynauds_SAS_Jacobs_24_11_2022
````

###############################
# Step 3 - liftover           #
###############################

````unix 

cd /genesandhealth/red/BenJacobs/raynaud_gwas/
Rscript ./scripts/liftover.R \
-g ./outputs/genes_and_health_raynauds_SAS_Jacobs_24_11_2022_raynaud.regenie \
-i hg38 \
-o hg19

cp ./liftover_temp_2022-11-24_27/gwas_results_hg19 ./outputs/

````
