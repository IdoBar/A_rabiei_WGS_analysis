# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# load packages we're going to use
pacman::p_load(tidyverse, vcfR, R.utils, RColorBrewer, ggrepel, glue)

# maybe reset genotype calls in which the Alternative allele depth is more than 5% than the Reference allele depth

#### Estimate error rates ####
analysis_name <- "Snippy_multi_07_07_2019"
analysis_folder <- file.path("..", analysis_name)
analysis_outdir <- "./output/duplicate_analysis/"
# load sequencing table
sequencing_table <- readxl::read_excel("./sample_info/A_rabiei_isolate_list_for_wgs.xlsx",
                                       sheet = "submission_info") %>% #arrange(Isolate) %>% 
  group_by(Isolate) %>% 
  mutate(counter=row_number(), Isolate_head=paste(Isolate, counter, sep="."))

snippy_sample_tbl <- read_tsv(file.path(analysis_folder, "snippy_input_samples.tab"),
                              col_names = c("snippy_name", "Read1", "Read2")) %>% 
  mutate(Submission_id=sub("trimmed_(.+)_R1.fastq.gz", "\\1", Read1))
snippy_dict <- setNames(snippy_sample_tbl$snippy_name, 
                         snippy_sample_tbl$Submission_id)
samples_dict <- setNames(sequencing_table$Isolate_head, sequencing_table$Submission_id)
source("./src/estimate_error_rates_vcf_files.R")
# load VCF file (Bowtie2 pipeline)
filtered_vcf_files <- file.path(analysis_folder, "core.vcf")
# error_rates <- tibble()
# for (f in filtered_vcf_files){
#   # f <- filtered_vcf_files[1]
#   # check error rates between replicate samples
#   
#   filtered_vcf <- read.vcfR(f)
#   vcf_headers <- colnames(filtered_vcf@gt)[-1]
#   sample_isolates <- samples_dict[vcf_headers]
#   colnames(filtered_vcf@gt)[-1] <- sample_isolates
#   fixed_vcf_file <- paste0(f, ".fixed.gz")
#   write.vcf(x=filtered_vcf,file = fixed_vcf_file)
#   R.utils::gunzip(fixed_vcf_file, overwrite=TRUE)
#   
#   error_rates <- rbind(error_rates, estimate_error_rates(paste0(f, ".fixed"), grouping_suffix = ".[0-9]+$"))
# }
error_rates <- estimate_error_rates(filtered_vcf_files[1], 
                                    grouping_suffix = "[abcdefg]") %>% 
write_xlsx(., excel_file = "output/results/pipeline_comparison.xlsx",
                     sheet=analysis_name, overwritefile=FALSE, asTable = TRUE )

# Check how to compare the error rates for each group.
# error_rates_sum <- error_rates %>% arrange(replicate_group) %>% mutate(analysis=str_extract(vcf_file, "FB_vars.+_2019")) %>% 
#   # mutate(analysis=sub("01_03", "BBtools_01_03", sub("07_04", "BT2_07_04", analysis))) %>%
#   dplyr::select(analysis, replicate_group, allele_error_rate) %>% 
#   # group_by(replicate_group) %>% 
#   spread(key=analysis, value = allele_error_rate) %>% 
#   mutate(allele_error_rate_diff=FB_vars_BBtools_01_03_2019-FB_vars_BT2_07_04_2019)
# openxlsx::write.xlsx(error_rates_sum, file = "output/results/pipeline_comparison.xlsx", asTable = TRUE)
# mean(error_rates_sum$allele_error_rate_diff)
# mean(error_rates_sum$FB_vars_07_04_2019)
# 
# #### Filtration ####
# analysis_basename="A_rabiei_2018"
# variant_method <- "haplo.bt2.fb"
# analysis_folder <- "../FB_vars_07_04_2019"

# f <- filtered_vcf_files[1]
# vcf <- read.vcfR(f)
# 
# # Check allelic depth (minor/major<0.15)
# 
# ad <- extract.gt(vcf, element = "AD", as.numeric = FALSE)
# check_ad <- function(allelic_depth, maf=0.15){
#   # allelic_depth=ad[1,9]
#   if (is.na(allelic_depth)) return(NA)
#   alleles <- as.numeric(unlist(str_split(allelic_depth, pattern = ",")))
#   minor_allele <- min(alleles)
#   major_allele <- max(alleles)
#   return(minor_allele/major_allele<maf)
# }
# 
# # check maf in markers
# test <- apply(ad, 2, function(a) a %>% purrr::map_lgl(check_ad))
# # check how many got reset from each genotype
# apply(test, 2, function(t) sum(!t, na.rm = TRUE))
# # clear all calls where the maf is too high
# vcf@gt[,-1][ test == FALSE ] <- NA
# 
# # clear depth that deviates from the quantiles (per sample)
# dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
# # apply(dp, 2, function(d) max(d, na.rm = TRUE))
# # find outliers
# min_depth <- 10
# max_depth <- 500
# # sample_outliers <- apply(dp, 2, function(d) sum(d>max_depth, na.rm = TRUE))
# # find markers with extremely high depth
# marker_outliers <- apply(dp, 1, function(d) sum(d>max_depth, na.rm = TRUE))
# marker_outliers[marker_outliers>1] 
# # find markers with extremely low depth
# marker_outliers <- apply(dp, 1, function(d) sum(d<min_depth, na.rm = TRUE))
# low_dp_markers <- marker_outliers[marker_outliers>20]
# # remove low-depth markers
# vcf <- vcf[!marker_outliers>20,]
# 
# # Remove missing markers with vcfR
# # Check missingness in genotypes
# # dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
# #
# 
# # Check missingness in markers
# markerMiss <- apply(gt, MARGIN = 1, function(x){ sum(is.na(x)) })
# markerMissRate <- markerMiss/ncol(vcf@gt[,-1])
# miss_thresh <- 0.25
# marker_missing <- markerMissRate[markerMissRate>miss_thresh]
# # remove markers with more than 25% missing 
# clean_vcf <- vcf[markerMissRate<miss_thresh,]
# hist(markerMissRate, col = "#8DD3C7", xlab = "Missingness (%)")
# gt <- extract.gt(clean_vcf)
# myMiss <- apply(gt, 2, function(x) sum(is.na(x)) )
# percMiss <- myMiss*100/nrow(clean_vcf)
# # plot missingness
# palette(brewer.pal(n=12, name = 'Set3'))
# # pdf(filedate(glue::glue("{analysis_basename}_missingness"), 
# #              ext = ".pdf", 
# #              outdir = glue("{analysis_outdir}/plots/")),
# #     width = 12, height = 9)
# # par(mar = c(12,4,4,2))
# barplot(percMiss, las = 2, col = 1:12)
# title(ylab = "Missingness (%)")
# # dev.off()
# # write back to file
# vcf_basename <- tools::file_path_sans_ext(f)
# write.vcf(clean_vcf, file = glue::glue("{vcf_basename}.clean.vcf.gz"))
# 

#### EDA ####
# filtered_vcf <- read.vcfR(filtered_vcf_file)
genind_obj <- vcfR2genind(filtered_vcf)
# Convert to allele frequencies (manage missing data)
X <- adegenet::tab(genind_obj, freq = TRUE, NA.method = "mean")
pca1 <- ade4::dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 20)

pca_data <- pca1$li %>% rownames_to_column("snippy_name")  %>%
  inner_join(snippy_sample_tbl %>% select(snippy_name, Submission_id)) %>%
  # mutate(seq_fontface=if_else(Sequenced=="Illumina", "plain", "bold"),
  #        Pathotype=factor(Pathotype, levels=path_levels),
  #        Coverage=round(Coverage, 2)) %>% 
  left_join(sequencing_table %>%  dplyr::select(-one_of(c("counter", "Comments"))))# %>% 
  # mutate(Sequencing_Centre=if_else(is.na(Sequencing_Centre), "AG", Sequencing_Centre))

# Calculate variance for each component (columns starting with 'C')
pca_var <- pca1$li %>% summarise_at(vars(starts_with("Axis")), var)
# Calculate the percentage of each component of the total variance
percentVar <- pca_var/sum(pca_var)
# Define colour palette, but get rid of the awful yellow - number 6
pal <- brewer.pal(9, "Set1")
pal2 <- brewer.pal(8, "Dark2")

seq_cols <- adjustcolor( c("grey15", "dodgerblue3"), alpha.f = 0.8)
# seq_face <- setNames(c("bold", "plain"),unique(pca_data$Sequenced))

shapes <- c(21:25)
# Isolate to annotate
replicated_samples <- pca_data %>% count(Isolate) %>% filter(n>1)
# outliers <- pca_data %>% filter(Axis2< -2.5) %>% .[,"Isolate"]
outliers <- pca_data %>% inner_join(replicated_samples) %>% .[,"snippy_name"]
# outliers <- c("16CUR017", "TR9543", "TR9529", "15CUR003") # PacBio-sequenced
# outliers <- pca_data$Isolate
#### PCA plot 1-2 biological ####
# shape - Year, Seq.Provider, Host (State)
# Colors - Haplotype (diverged), State (diverged), Pathotype (RdYlGn)
# Size - Pathogenicity
# Create the plot for C1 and C2
p <- ggplot(pca_data, aes(x=Axis3, y=Axis4, color=Sequencing_Centre, 
                          # size=Year,
                          shape=Sequencing_Centre)) # Host
# Plot and add labels (represent component variance as a percentage)
p + geom_point(alpha = 0.85, size=5) + # , stroke = 1
  # geom_text_repel(aes(colour=Sequenced, label=Isolate, fontface=seq_fontface), show.legend = FALSE,size=3.5, point.padding = 0.75) +
  # geom_text_repel(aes(label=Isolate_head, colour=adjustcolor( c("dodgerblue3"), 
  #                                                             alpha.f = 0.8)),
  #             show.legend = FALSE, size=3.5, point.padding = 0.75) +
  # scale_color_() +
  scale_color_brewer(palette = "Set1", direction = 1) +
  # scale_fill_manual(values = pal, aesthetics = c("fill")) + 
  # scale_shape_manual(values = shapes) +
  # scale_color_manual(values=seq_cols) +
  # scale_size_continuous(range = c(3,6)) +
  guides(shape = guide_legend(override.aes = list(size = 5 , color=pal[1:3]), 
                              order = 1, title = "Seq.centre"), #
         color = "none") +   #  guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2)
  # scale_size_manual(values = path_sizes[pca_data$Pathogenicity]) +
  
  # geom_text_repel(aes(label=Haplotype), size=4, point.padding = 0.75) +
  plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue::glue("C3: {round(percentVar[3]*100,2)}% variance"),
        y=glue::glue("C4: {round(percentVar[4]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("{analysis_name}_Ind_PCA_seq_centre_PC3-4"),
                ext = ".pdf", outdir = glue::glue("{analysis_outdir}/plots"),
                dateformat = FALSE), width = 10, height=8)
save.image(filedate(glue("{analysis_name}_duplicate_analysis"), ext = ".RData", glue("{analysis_outdir}")))

