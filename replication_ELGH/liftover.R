# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
options(warn = -1)

# args
option_list = list(
  make_option(c("-g","--gwas_path"), type = "character", default = "in", help = "path to gwas file", metavar = "gwas"),
  make_option(c("-i","--input_build"), type = "character", default = "hg38", help = "input build", metavar = "input"),
  make_option(c("-o","--output_build"), type = "character", default = "hg19", help = "output build", metavar = "input")
)

opt = parse_args(OptionParser(option_list = option_list))

# arg checking
if(!file.exists(opt$gwas_path)){
  stop(opt$gwas_path, " must be a valid path to REGENIE GWAS results")
}
if(!(opt$input_build %in% c("hg19","hg38")) | !(opt$output_build %in% c("hg19","hg38")) ){
  stop("Both input and output build must be one of hg19 or hg38.")
}

gwas_path = opt$gwas_path

# read in chainfile

message("Converting from ",opt$input_build," to ",opt$output_build)
chainfile_path = if(opt$input_build == "hg38" & opt$output_build =="hg19"){
  "/genesandhealth/library-green/sanger/Liftover_chain/hg38ToHg19.over.chain.gz"
} else if(opt$input_build == "hg19" & opt$output_build =="hg38"){
  "/genesandhealth/library-green/sanger/Liftover_chain/hg19ToHg38.over.chain.gz"
}


# read in gwas
message("Reading in GWAS")
gwas_res = read_table(file = gwas_path, col_types = cols(.default = "c"))

# convert to bed
message("Converting to bed file format")
gwas_bed = gwas_res %>%
  mutate(newchr = paste0("chr",CHROM),
         newpos = as.numeric(GENPOS) - 1,
         newpos2 = GENPOS,
         originalid = ID) %>%
  dplyr::select(newchr,newpos,newpos2,originalid)

# write to file
today = Sys.Date()
rnum = as.integer(runif(min=1,max=1000,n=1))
dir_to_make = paste0("liftover_temp_",today,"_",rnum)
message("Results being written to temporary folder: ",dir_to_make)
system(paste0("mkdir ",dir_to_make))
write_tsv(gwas_bed,paste0(dir_to_make,"/gwas_bed.bed"),col_names = F)

# run liftover
cmd = paste0("liftOver ",dir_to_make,"/gwas_bed.bed ",chainfile_path," ",dir_to_make,"/gwas_output ",dir_to_make,"/unlifted")
system(cmd)

# get lifted coords
new_coords = read_table(paste0(dir_to_make,"/gwas_output"),
                        col_names = F,
                        cols(.default = "c"))
colnames(new_coords) = c("CHR","BP_1","BP","ID")
new_coords = new_coords %>%
  dplyr::select(ID,BP) %>%
  distinct(ID,.keep_all = T)

# join
# remove duplicated positions
gwas_res = gwas_res %>%
  distinct(ID,.keep_all = T)

message(nrow(gwas_res), " SNPs in input GWAS")
gwas_res_lifted = gwas_res %>%
  filter(ID %in% new_coords$ID) %>%
  left_join(new_coords,by="ID") %>%
  mutate(old_GENPOS = GENPOS, old_ID = ID) %>%
  mutate(GENPOS = BP) %>%
  mutate(ID = paste0(CHROM,":",GENPOS)) %>%
  dplyr::select(-BP) %>%
  filter(!is.na(ID)) %>%
  dplyr::select(-old_ID,-old_GENPOS)
message(nrow(gwas_res_lifted), "SNPs successfully lifted to ",opt$output_build)

# write to file
outpath = paste0(dir_to_make,"/gwas_results_",opt$output_build)
write_tsv(gwas_res_lifted,file=outpath)
message("Writing output to: ",outpath)
message("All done")
