# devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# pacman::p_load(tidyverse, RColorBrewer, doParallel)
#
# #### Process vcf files ####
# analysis_basename="A_rabiei_2018_35_isolates"
# vcf_file <- recent_file("./output/", glue::glue("{analysis_basename}.+.vcf"))
#
# vcf_files <- dir(pattern = glue::glue("{analysis_basename}.+.vcf"), recursive = TRUE)
# # vcf_files <- vcf_files[grepl("stacks_.+/populations.r80/", vcf_files)]
# registerDoParallel(cores = 3)
#
# # Load vcf file
# vcf_file <- "./output/A_rabiei_2018_35_isolates_DP_corr.recal.filt.snps.vcf"
estimate_error_rates <- function(vcf_input, grouping_suffix){
  error_rates_df <- tibble(vcf_file=character(),
                               replicate_group=character(),
                               locus_error_rate=numeric(),
                               allele_error_rate=numeric())
                               # stringsAsFactors=FALSE)
  # i=1
  vcf_data <- read_tsv(vcf_input, comment = "##") #%>%
  # stacks_params <- stringr::str_match(vcf_files[i], "stacks_M(\\d)m(\\d)n(\\d)")
  # mutate(`#CHROM`=as.numeric(unlist(strsplit(`#CHROM`, split="_"))[2]))  # n_max=100
  sample_cols <- colnames(vcf_data)[10:ncol(vcf_data)]
  # Identify replicates
  sample_grouping <- tibble(sample=sample_cols,
                            replicate_group=gsub(grouping_suffix, "", sample_cols))
  groups <- sample_grouping %>% count(replicate_group) %>% filter(n>1) %>%
    .$replicate_group
  for (r in groups) {
    # r=groups[1]
    LogMsg(sprintf("Processing replicate group %s of file %s", r, vcf_input))
    replicates <- sample_grouping %>%
      filter(replicate_group==r)

    marker_errors <- apply(vcf_data[replicates$sample], 1, function(g) {
      # res_df <- data.frame(dropped_locus=FALSE,
      #                      missing_locus=FALSE,
      #                      allele_mismatch=FALSE,
      #                      stringsAsFactors=FALSE)
      locus_errors <- setNames(rep(FALSE, 3), c("dropped_locus", "missing_locus",
                                                "allele_mismatch"))
      # g <- vcf_data[7, replicates$sample]
      genotypes <- sub(":.+", "", g) %>% sub("0/1", "1/0", ., fixed = TRUE)
      # genotypes <- sub("1/0", "0/1", g)
      # log_vec <- ifelse(biallelic, length(prop_names[!grepl("\\.", prop_names)])==2,
      # length(prop_names[!grepl("\\.", prop_names)])>1)
      p_table <- prop.table(table(genotypes))
      # g_table <- p_table[!grepl("\\./\\.", names(p_table))]
      if (any(grepl("\\./\\.", names(p_table)))) {
        if (length(p_table)==1) {
          # res_df$dropped_locus[1] <- TRUE
          locus_errors["dropped_locus"] <- TRUE
        } else locus_errors["missing_locus"] <- TRUE # res_df$missing_locus[1] <- TRUE
        # res <- g_table[1]>0.05 && g_table[2]>0.05
      }  else {
        if (length(p_table)>1) locus_errors["allele_mismatch"] <- TRUE # res_df$allele_mismatch[1] <- TRUE
      }
      return(matrix(locus_errors, nrow = 1))
    }) %>% t() %>% as_tibble()
    colnames(marker_errors) <- c("dropped_locus", "missing_locus",
                                 "allele_mismatch")
    # marker_errors[1:10, 1:3]
    locus_error_rate <- sum(marker_errors$missing_locus)/nrow(marker_errors)
    allele_error_rate <- sum(marker_errors$allele_mismatch)/(nrow(marker_errors)-sum(marker_errors$dropped_locus))

    error_rates_df <- tibble("vcf_file"=vcf_input, "replicate_group"=r, "replicate_num"=nrow(replicates),
                                 "locus_error_rate"=locus_error_rate,
                                 "allele_error_rate"=allele_error_rate) %>%
      bind_rows(error_rates_df, .)

  }
  return(error_rates_df)
}

# head(error_rates_results)
# # str(error_rates_results)
# write_tsv(error_rates_results, filedate("stacks_param_error_rates", ".tsv", "data", FALSE))
