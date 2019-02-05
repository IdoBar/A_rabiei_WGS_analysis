vcf_filtration <- function(vcf_file, miss_rates = seq(0.2, 0.1, -0.05), geno_miss_rate=0.2, check_poly_parents="", remove_hetero=FALSE, remove_multi_allele=TRUE, write_filtered_vcf="", poly_only=TRUE) {
  LogMsg(sprintf("Filtering vcf file (%s)", vcf_file))
  vcf_data <- readr::read_tsv(vcf_file, comment = "##") #%>%
  # mutate(`#CHROM`=as.numeric(unlist(strsplit(`#CHROM`, split="_"))[2]))  # n_max=100
  sample_cols <- colnames(vcf_data)[10:ncol(vcf_data)]

  if (isTRUE(remove_multi_allele)) {
    multi_allelic_snps <- vcf_data %>% dplyr::filter(nchar(ALT)>1)
    LogMsg(sprintf("Filtered %d multi-allelic markers", nrow(multi_allelic_snps)))
    vcf_data <- vcf_data %>% dplyr::filter(nchar(ALT)==1)
  }

  if (isTRUE(remove_hetero)){
    # Relace missing genotypes and heterozygotes with ./.
    vcf_data <- vcf_data %>% dplyr::mutate_at(vars(one_of(sample_cols)),
                  .funs = funs(if_else(grepl("^1/0", .) | grepl("^0/1", .), "./.", .)))
    LogMsg(sprintf("Resetting heterozygote calls (to ./.)"))
  }
  vcf_data <- vcf_data %>% dplyr::mutate_at(vars(one_of(sample_cols)),
                                    .funs = funs(if_else(grepl("^\\.$", .), "./.", .)))
  marker_msg <- "markers"
  if (isTRUE(poly_only)){
    parents <- sample_cols[grepl(check_poly_parents, sample_cols)]
    poly_markers <- apply(vcf_data[parents], 1, function(g) {
      # g <- vcf_data[1, parents]
      genotypes <- sub(":.+", "", g) %>% sub("0/1", "1/0", ., fixed = TRUE)
      # log_vec <- ifelse(biallelic, length(prop_names[!grepl("\\.", prop_names)])==2,
      # length(prop_names[!grepl("\\.", prop_names)])>1)
      p_table <- prop.table(table(genotypes))
      g_table <- p_table[!grepl("\\./\\.", names(p_table))]
      if (length(g_table)>=2) {
        res <- TRUE # g_table[1]>0.05 && g_table[2]>0.05
      }  else res <- FALSE
      return(res)
    })

    vcf_data <- vcf_data %>% dplyr::filter(poly_markers)
    LogMsg(sprintf("Filtered %d non-polymorphic markers", sum(!poly_markers)))
    marker_msg <- paste("polymorphic", marker_msg)
  }
  # Check if we have any heterozygote left
  # sum(apply(vcf_data[sample_cols], 1, function(r) sum(grepl("^1/0", r) || grepl("^0/1", r))))
  exclude_markers <- FALSE
  vcf_filt_data <- vcf_data
  incl_samples <- sample_cols
  counter=1

  for (m in miss_rates){

    # sample_col_filt <- colnames(vcf_filt_data)[10:ncol(vcf_filt_data)]
    # m <- miss_rates[3]

    excl_samples <- incl_samples

    while (length(excl_samples)>0) {
      # vcf_filt_data <- vcf_filt_data %>% filter(!exclude_markers)
      exclude_markers <- apply(vcf_filt_data[incl_samples], 1, function(marker) {
        sum(grepl("^\\./\\.", marker) | grepl("^\\.$", marker)) > m*length(incl_samples)
      } )
      table(exclude_markers)

      vcf_filt_data <- vcf_filt_data %>% dplyr::filter(!exclude_markers)
      # Find samples with too many missing markers

      excl_samples <- incl_samples[apply(vcf_filt_data[incl_samples], 2,
                         function(geno) sum(grepl("\\./\\.", geno))>(geno_miss_rate)*nrow(vcf_filt_data))]


      # Now repeat marker exclusion without the problematic genotypes
      incl_samples <- incl_samples[!incl_samples %in% excl_samples]
      LogMsg(sprintf("Genotype and Marker filtering, round %d, missing rate %.2f", counter, m))
      LogMsg(sprintf("Filtered %d genotypes and %d markers", length(excl_samples), sum(exclude_markers)))
      counter <- counter + 1
    }

  }

  exclude_genotypes <- sample_cols[!sample_cols %in% incl_samples]
  # Remove missing genotypes and markers, update marker information (CHROM, POS, LOC)
  clean_vcf <- vcf_filt_data %>%
    dplyr::select(-one_of(exclude_genotypes))

  # table(poly_markers)

  LogMsg(sprintf("Filtering VCF completed: %d %s from %d genotypes remained", nrow(clean_vcf), marker_msg, ncol(clean_vcf)-9))

  # vcf_filt_data <- vcf_filt_data %>% filter(poly_markers)
   #%>%
  # Write back vcf
  vcf_basename <- tools::file_path_sans_ext(basename(vcf_file))
  if (write_filtered_vcf=="") return(clean_vcf)
  clean_vcf_file <- ifelse(write_filtered_vcf=="default", filedate(sprintf("%s.Marker_miss%.0fGeno_miss%.0f",vcf_basename, m*100,geno_miss_rate*100), ".clean.vcf","output", dateformat = FALSE), write_filtered_vcf)


  # sub(".clean.sorted.vcf", "_miss40.clean.vcf",vcf_file, fixed = TRUE))

  # unlink(clean_vcf_file)
  # Copy the header into a new file
  system2("grep", args=c("'^##'", vcf_file), stdout = clean_vcf_file)
  if (file.exists(clean_vcf_file)) readr::write_tsv(clean_vcf, clean_vcf_file, append = TRUE, col_names = TRUE)

  # completed_tasks <<- 1
  return(clean_vcf)
}
