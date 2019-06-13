# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# load packages we're going to use
pacman::p_load(tidyverse, vcfR, R.utils)

# maybe reset genotype calls in which the Alternative allele depth is more than 5% than the Reference allele depth

#### Estimate error rates ####
# load sequencing table
sequencing_table <- readxl::read_excel("./sample_info/A_rabiei_isolate_list_for_wgs.xlsx",
                                       sheet = "submission_info") %>% #arrange(Isolate) %>% 
  group_by(Isolate) %>% 
  mutate(counter=row_number(), Isolate_head=paste(Isolate, counter, sep="."))

samples_dict <- setNames(sequencing_table$Isolate_head, sequencing_table$Submission_id)
source("./src/estimate_error_rates_vcf_files.R")
# load VCF file (Bowtie2 pipeline)
filtered_vcf_files <- c("../FB_vars_07_04_2019/A_rabiei_2018_isolates_haplo_ind.bt2.fb.Qual10.U20000DP.poly.snps.vcf", "../FB_vars_01_03_2019/A_rabiei_2018_isolates_haplo_ind.bbtools.fb.Qual10.U20000DP.poly.snps.vcf")
error_rates <- tibble()
for (f in filtered_vcf_files){
  # check error rates between replicate samples
  
  filtered_vcf <- read.vcfR(f)
  vcf_headers <- colnames(filtered_vcf@gt)[-1]
  sample_isolates <- samples_dict[vcf_headers]
  colnames(filtered_vcf@gt)[-1] <- sample_isolates
  fixed_vcf_file <- paste0(f, ".fixed.gz")
  write.vcf(x=filtered_vcf,file = fixed_vcf_file)
  R.utils::gunzip(fixed_vcf_file, overwrite=TRUE)
  
  error_rates <- rbind(error_rates, estimate_error_rates(paste0(f, ".fixed"), grouping_suffix = ".[0-9]+$"))
}
# Check how to compare the error rates for each group.
error_rates_sum <- error_rates %>% arrange(replicate_group) %>% mutate(analysis=str_extract(vcf_file, "FB_vars.+_2019")) %>% 
  mutate(analysis=sub("01_03", "BBtools_01_03", sub("07_04", "BT2_07_04", analysis))) %>% 
  dplyr::select(analysis, replicate_group, allele_error_rate) %>% 
  # group_by(replicate_group) %>% 
  spread(key=analysis, value = allele_error_rate) %>% 
  mutate(allele_error_rate_diff=FB_vars_BBtools_01_03_2019-FB_vars_BT2_07_04_2019) # %>% 
 
openxlsx::write.xlsx(error_rates_sum, file = "output/results/pipeline_comparison_test.xlsx", asTable = TRUE)

mean(error_rates_sum$allele_error_rate_diff)
mean(error_rates_sum$FB_vars_07_04_2019)
