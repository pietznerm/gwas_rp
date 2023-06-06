#!/bin/sh
## Author Dr. Summaira Yasmeen

## script to run the GWAS 
## Requires regenie_v2.2.4
## HPC SLURM

################ Specify SLURM JOB needs
#SBATCH --job-name GWAS_raynaud
#SBATCH --partition=compute
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 120
#SBATCH --array=1-24%5
#SBATCH --mail-type=FAIL
#SBATCH --output=~../log/REGENIE_GWAS-%j.log

######################Export relevant directories 
## general 
export ukbb_data=.../data/UK_biobank
## as data directories
export geno=${ukbb_data}/genotypes/raw_data
export imputed_data=${ukbb_data}/genotypes/bgen_files
export imputed_qc_snps=${ukbb_data}/genotypes/variant_qc

## anaylsis directory..
export current_proj=.../Raynaud_GWAS
export input=${current_proj}/input
export tmpdir=${current_proj}/tmpdir
export output=${current_proj}/output

##### Create variable names 

 ## SNPs that passed QC
 snplist_step1="qcPASS.snplist_raw_genotypes.txt"
 ## IDS that pass QC
 gen_ids="UKBB_array_data_individuals_with_genotypes.txt"
 ## Covariates
 cov_file="covariates.txt"
 ## Phenotype
 pheno="phenotype.txt" 

 ## Variable to store output from step 1
 step1_output="ukb.${pheno}_step1"
 tmpfile="regenie_tmp_preds_${pheno}"
 ## Samples at Chr X
 sample_x="ukb_imp_chrX.sample"
 ## Sample at all other Chromsome
 sample_all="ukb_imp_chr1.sample"
 
 ## Variable to store output at step2
 step2_output="step1_chr${chr}_"

## get the chromosome
echo "Job ID: $SLURM_ARRAY_TASK_ID"


## STEP1
/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
--step 1 \
--bed ${geno}/ukb_cal_allChrs \
--extract ${input}/${snplist_step1} \
--keep ${input}/${gen_ids} \
--phenoFile ${input}/${pheno} \
--covarFile ${input}/${cov_file} \
--threads 30 \
--bt \
--bsize 1000 \
--lowmem \
--lowmem-prefix ${tmpdir}/${tmpfile} \
--out ${input}/${step1_output}


## Perform Step 2

## get the chromosome
if [[ $chr -eq 23 ]]; then
 export chr="X"
echo "echo $pop.. $chr .. $batch"

## create variant inclusion list for the respective chromosome

 cat ${var}/ukb_imp_chr${chr}_qced.txt | awk -v chr=${chr} '{if(NR != 1 && $3 == chr) print $2}' - > ${tmpdir}/tmp.ex.chr${chr}.list

/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
	--step 2 \
 	--bgen ${imputed_data}/ukb_imp_chr${chr}.bgen \
 	--ref-first \
 	--extract ${tmpdir}/tmp.ex.chr${chr}.list \
 	--sample  ${imputed_data}/${samplex} \
 	--phenoFile ${input}/phen/${pheno} \
 	--covarFile ${input}/cov/${cov_file} \
 	--threads 30 \
 	--bt \
 	--spa 0.01 \
 	--pred ${input}/${step1_output}_pred_list \
 	--bsize 400 \
 	--out ${output}/${step2_file}

## remove variant inclusion list by chr number..
rm ${tmpdir}/tmp.${pop}.ex.chr${chr}.list

else
echo "echo  $chr .. "
  
  ## create variant inclusion list for the respective chromosome
 cat ${var}/ukb_imp_${pop}_chr${chr}_qced.txt | awk -v chr=${chr} '{if(NR != 1 && $3 == chr) print $2}' - > ${tmpdir}/tmp.ex.chr${chr}.list

/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
	--step 2 \
	--bgen ${imputed_data}/ukb_imp_chr${chr}.bgen \
	--ref-first \
	--extract ${tmpdir}/tmp.ex.chr${chr}.list \
	--sample  ${imputed_data}/${sample_all} \
	--phenoFile ${input}/${pheno} \
	--covarFile ${input}/${cov_file} \
	--threads 30 \
	--bt \
	--spa 0.01 \
	--pred ${input}/${step1_output}_pred_list \
	--bsize 400 \
	--out ${output}/${step2_file}

## remove variant inclusion list
rm ${tmpdir}/tmp.ex.chr${chr}.list
fi
## change to relevant directory
