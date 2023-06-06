#!/bin/sh
## Author Dr. Summaira Yasmeen

## script to run compute genetic correlations
## Requires LDSC: https://github.com/bulik/ldsc
## HPC SLURM
## script to runs Genetic Correlation between Raynaud GWAS and phecodes + selected traits

#SBATCH --job-name Genetic Correlation
#SBATCH --partition=compute
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 120
#SBATCH --array=1-30%5
#SBATCH --mail-type=FAIL
#SBATCH --output=~../log/Genetic Correlation-%j.log



## project in running dir..
export current_proj=.../genetic_correlations
export input=${current_proj}/input
export output=${current_proj}/output

## output 
export tmpdir=${current_proj}/tmpdir
export LDSC=../ldsc


selected_traits="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' selected_traits.list.txt)"
## specify the rg string 
## example raynaud_gwas,phe1_gwas,phe2_gwas,...,pheN_gwas
rg_string="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' rg.string)"

echo "Perform Gene_Correlation for ${selected_traits} .. Job ID: $SLURM_ARRAY_TASK_ID"

## -->  Munge summary stats <-- ##
${LDSC}/munge_sumstats.py \
	--sumstats ${tmpdir}/${selected_traits} \
	--snp SNP \
	--a1 REF \
	--a2 ALT \
	--p Pval \
	--signed-sumstats BETA,0 \
	--out ${tmpdir}/${selected_traits} \
	--chunksize 500000 \
	--keep-maf \
	--merge-alleles ${LDSC}/w_hm3.snplist

## --> run gene correlation<-- ##

/ldsc.py \
	--rg  ${rg_string} \
	--ref-ld-chr ${LDSC}/eur_w_ld_chr/ \
	--w-ld-chr ${LDSC}/eur_w_ld_chr/ \
	--out ${output}/${selected_traits}_rg



