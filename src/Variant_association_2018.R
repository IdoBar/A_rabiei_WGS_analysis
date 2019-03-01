devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
# Install these to enable installation of radiator and grur
# install_cran(c("future", "Rcpp", "apex", "copula", "gsl", "ADGofTest", "stabledist", "pcaPP", "pspline", "shinyFiles",
#                "RJSONIO", "swfscMisc", "mapdata", "maps"))
# install.deps(c("tidyverse/glue", "tidyverse/tidyselect", "thierrygosselin/radiator", "thierrygosselin/grur"), repo="git") # need to be run first, otherwise fails to re-install tidyverse packages if they're already loaded
# devtools::install_github("thierrygosselin/radiator")
# install.deps(c("thierrygosselin/radiator"), repo="git")
# Install CRAN archived packages
devtools::install_version("GenABEL.data",version="1.0.0")
devtools::install_version("GenABEL",version="1.8-0")
CRAN_packages <- c("tidyverse", "RColorBrewer", "ggrepel","GenABEL", "outliers", "adegenet", "SNPassoc", "vcfR", "glue")
pacman::p_load(char=CRAN_packages)
# install.packages("SNPassoc")
# install.deps(CRAN_packages)

# install.deps("SNPassoc", repo = "bioc")
# Download working version of SNPassoc from http://www.creal.cat/media/upload/arxius/jr/SNPassoc/SNPassoc_1.8-5.zip and copy to R library

#### Read data files ####
# Read sample metadata and save to file
analysis_basename="A_rabiei_2018"
analysis_folder <- "../BB_vars_26_02_2019"
sequencing_table <- readxl::read_excel("./sample_info/A_rabiei_isolate_list_for_wgs.xlsx",sheet = "submission_info")
sequencing_dict <- set_names(sequencing_table$Isolate, sequencing_table$Submission_id)
samples_table <- readxl::read_excel("./sample_info/A_rabiei_isolate_list_for_wgs.xlsx") #,
                                    # sheet = "Sequenced")

# Read mapping stats (after preparation by the markdown notebook)
mapping_stats <- recent_file("./output/results", glue("{analysis_basename}.+mapping.stats.txt")) %>%
  read_tsv() 
mapping_summary <- recent_file("./output/results", glue("{analysis_basename}.+mapping.sum.txt")) %>%
  read_tsv() 
# mapping_summary <- mapping_stats %>% filter(Coverage>20) %>%
#   group_by(Sequencing_Centre) %>%
#   summarise(Coverage=sprintf("x%.2f",mean(Coverage)),
#             Mapping_rate=sprintf("%.2f%%", mean(Mapping_rate, na.rm=TRUE)*100),
#             Mapping_qual=mean(Mapping_quality_mean),
#             Files=n()) %>%
#   mutate(Reads_per_file=format(c(36500108.00, 2224642.64),
#                                digits = 2, scientific = TRUE)) %>%
#   arrange(desc(Sequencing_Centre))

# Save as a tfam file

path_levels <- c("Low", "Medium", "Moderate", "High", "Very High")
strata <- samples_table %>% # filter(!is.na(Sequencing_Centre)) %>%
  mutate(INDIVIDUALS=Isolate, STRATA=path_levels[Pathogenicity+1],
         Site=paste(Site, State, sep=", ")) %>%
  dplyr::select(INDIVIDUALS, STRATA, Pathogenicity, Site, State, Collection_Year) %>%
  left_join(mapping_stats[c("Isolate", "Coverage")],
                                      by= c("INDIVIDUALS"="Isolate")) # %>% write_tsv(filedate(analysis_basename, ".strata", "./output/results", dateformat = FALSE))


#### VCF EDA with vcfR ####
# vcf_file <- recent_file("./data", glue::glue("{analysis_basename}.+.vcf"))
vcf_file <- recent_file(analysis_folder, glue("{analysis_basename}_isolates.+.vcf"))
vcf_basename <- tools::file_path_sans_ext(basename(vcf_file))
vcf <- read.vcfR(vcf_file)


# Fix sample_names
# colnames(vcf@gt)[colnames(vcf@gt) %in% agrf_table$Submission_ID] <- agrf_dict[colnames(vcf@gt)[colnames(vcf@gt) %in% agrf_table$Submission_ID]]


#### QA and EDA ####
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
head(dp)

# plot depth
dpf <- reshape2::melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
dpf <- dpf[ dpf$Depth > 0,]
p <- ggplot(dpf, aes(x=Sample, y=Depth)) + geom_violin(fill="#C0C0C0", adjust=1.0,
                                                       scale = "count", trim=TRUE)
p + theme_bw() + theme(axis.title.x = element_blank(),
               axis.text.x = element_text(angle = 60, hjust = 1, size=12)) +
#  p <- p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
  scale_y_continuous(trans=scales::log2_trans(),
                            breaks=c(1, 10, 100, 800),
                            minor_breaks=c(1:10, 2:10*10, 2:8*100)) +
  theme(axis.title.y = element_text(size=12), panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6), panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2)) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2)
ggsave(filename = glue::glue("./output/plots/{analysis_basename}_depth_dist.pdf"),
       width = 12, height = 8)

# Check missingness within samples
head(extract.gt(vcf))
gt <- extract.gt(vcf)
myMiss <- apply(gt, 2, function(x){ sum(is.na(x)) })
percMiss <- myMiss*100/nrow(vcf)

palette(brewer.pal(n=12, name = 'Set3'))

par(mar = c(12,4,4,2))
barplot(percMiss, las = 2, col = 1:12)
title(ylab = "Missingness (%)")

# Check missingness in markers
markerMiss <- apply(gt, MARGIN = 1, function(x){ sum(is.na(x)) })
markerMissRate <- markerMiss/ncol(vcf@gt[,-1])

hist(markerMissRate, col = "#8DD3C7", xlab = "Missingness (%)")
# clear depth that deviates from the quantiles (per sample)
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
# Visualise depth
pdf(glue::glue("./output/plots/{analysis_basename}_depth_heatmap.pdf"),
    width = 12, height = 9)
heatmap.bp(dp, rlabels = FALSE)
dev.off()

#### Filtration ####
# find outliers
min_depth <- 5
max_depth <- 500
sample_outliers <- apply(dp, 2, function(d) sum(d>max_depth))
marker_outliers <- apply(dp, 1, function(d) sum(d>max_depth))

marker_outliers[marker_outliers>20]
max(marker_outliers)

sum(dp>max_depth)
sum(dp<min_depth)
# quants <- apply(dp, MARGIN=2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)
# dp2 <- sweep(dp, MARGIN=2, FUN = "-", quants[1,])
# dp[dp2 < 0] <- NA
# dp2 <- sweep(dp, MARGIN=2, FUN = "-", quants[2,])
# dp[dp2 > 0] <- NA
# clear depth by minimum
dp[dp < min_depth] <- NA
# clear depth by maximum
dp[dp > max_depth] <- NA

sum(is.na(dp))
# clear all calls where the depth is too low or high
vcf@gt[,-1][ is.na(dp) == TRUE ] <- NA



# Remove missing markers with vcfR
# Check missingness in markers
# dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
#
gt <- extract.gt(vcf)
myMiss <- apply(gt, 2, function(x){ sum(is.na(x)) })
percMiss <- myMiss*100/nrow(vcf)

palette(brewer.pal(n=12, name = 'Set3'))

par(mar = c(12,4,4,2))
barplot(percMiss, las = 2, col = 1:12)
title(ylab = "Missingness (%)")

# Check missingness in markers
markerMiss <- apply(gt, MARGIN = 1, function(x){ sum(is.na(x)) })
markerMissRate <- markerMiss/ncol(vcf@gt[,-1])

hist(markerMissRate, col = "#8DD3C7", xlab = "Missingness (%)")
write.vcf(vcf, file = glue::glue("./output/data_files/{vcf_basename}.DP_corr.vcf.gz"))


#
# filtered_vcf <- vcf[markerMissRate < 0.1, ]
# # Save filtered file
# write.vcf(filtered_vcf, "./output/A_rabiei_2018_35_isolates_DP_filt.recal.vcf")
#
# # Visualise filtered set
# dp2 <- extract.gt(filtered_vcf, element = "DP", as.numeric=TRUE)
# pdf(glue::glue("./plots/{analysis_basename}_miss_cov_1000_heatmap.pdf"),
#     width = 12, height = 9)
# heatmap.bp(dp2[1:1000,], rlabels = FALSE)
# dev.off()

# custom filtration
source("./src/vcf_filtration.R")
snpsift_filtered_file <- recent_file(analysis_folder, glue("{analysis_basename}_isolates.+passed.snps.vcf"))
vcf_basename <- tools::file_path_sans_ext(basename(snpsift_filtered_file))
filtered_vcf_file <-  glue::glue("./output/data_files/{vcf_basename}.poly.vcf")
vcf_tab_filt <- vcf_filtration(snpsift_filtered_file, miss_rates = seq(0.2,0.1, -0.05),
            geno_miss_rate=0.1, remove_hetero = FALSE, remove_multi_allele = FALSE, write_filtered_vcf = "", poly_only = TRUE)


clean_vcf_data <- read_tsv(filtered_vcf_file, comment = "##")
sample_cols <- colnames(clean_vcf_data)[10:ncol(clean_vcf_data)]
# Fix sample names
AGRF_samples <- which(colnames(vcf_tab_filt) %in% agrf_table$Submission_ID)
if (length(AGRF_samples)>0) colnames(vcf_tab_filt)[AGRF_samples] <- agrf_dict[colnames(vcf_tab_filt)[AGRF_samples]]

# output as tfam format
tfam_table <- samples_table %>% filter(!is.na(Sequencing_Centre),
                                       Isolate %in% sample_cols) %>%
  dplyr::select(State, Isolate) %>%
  mutate(V3=0, V4=0, V5=0, V6=-9) %>%
  write_tsv(filedate(analysis_basename, ".tfam", "./output", dateformat = FALSE), col_names = FALSE)

# write filtered vcf to file (not needed if done through the function)
# system2("grep", args=c("'^##'", vcf_file), stdout = filtered_vcf_file)
# if (file.exists(filtered_vcf_file)) readr::write_tsv(vcf_tab_filt, filtered_vcf_file, append = TRUE, col_names = TRUE)


# check error rates between replicate samples
source("./src/estimate_error_rates_vcf_files.R")
error_rates <- estimate_error_rates(filtered_vcf_file, grouping_suffix = "-[A-z]+$")
# read genome data
dna_file <- "../../A_rabiei_me14_short_names.fasta"
dna <- ape::read.dna(dna_file, format = "fasta")
gff_file <- "../../Arab_me14_short_names.gff"
gff <- read.table(gff_file, sep="\t", quote="")

# chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
# chrom <- proc.chromR(chrom, verbose = TRUE)
#
# head(chrom)
# plot(chrom)
# chromoqc(chrom, dp.alpha = 22)
# chrom <- masker(vcf, min_QUAL=20, min_DP=150, max_DP=2000)
# chrom <- proc.chromR(chrom, verbose = FALSE)
# plot(chrom)

#### Analysis ####
filtered_vcf <- read.vcfR(filtered_vcf_file)
genind_obj <- vcfR2genind(filtered_vcf)
# load to adegenet genlight (straight from tped file)
# obj <- genomic_converter("data/A_rabiei_on_me14_DP_GQ_corr.bt2.fb.vcf", output = "genind", strata = paste0(analysis_basename, ".strata"),  snp.ld = "maf",parallel.core = 0) #  imputation.method="rf",
# # fix chromosome names
# obj$tidier_data <- obj$tidy.data %>% mutate(MARKERS=sub("Arab_Me14_", "", MARKERS), CHROM=sub("Arab_Me14_", "", CHROM))

# Analyze missingness
# ibm.rabiei <- missing_visualization("data/A_rabiei_on_me14.bt2.fb.vcf", strata = "data/A_rabiei_isolates.strata", strata.select = c("STRATA", "Site"))
# genind_obj <- obj$genind.no.imputation

#### Perform some EDA  ####
# Convert to allele frequencies (manage missing data)
X <- tab(genind_obj, freq = TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 4)

pca_data <- pca1$li %>% rownames_to_column("Isolate")  %>%
  left_join(., strata, by = c( "Isolate"= "INDIVIDUALS")) %>%
  mutate(Pathotype=factor(STRATA, levels=path_levels),
      Sequencing_Platform=if_else(Sequencing_Centre=="AGRF", "NextSeq", "HiSeq2500")) %>%
  filter(!is.na(Pathotype))

# Calculate variance for each component (columns starting with 'C')
pca_var <- pca1$li %>% summarise_at(vars(starts_with("Axis")), var)
# Calculate the percentage of each component of the total variance
percentVar <- pca_var/sum(pca_var)
# Define colour palette, but get rid of the awful yellow - number 6
pal <- brewer.pal(9, "Set1")[-6]
# path_guide <- tibble(path_levels, col=pal[c(3,2,5,1,4)], size=seq(2,4,length.out = 5))
path_sizes <- setNames(seq(2,4,length.out = 5), path_levels)
path_cols <- setNames(pal[c(3,2,5,1,4)], path_levels)
# Isolate to annotate
outliers <- c("TR9571", "TR9543", "TR9529", "15CUR003")
# Create the plot for C1 and C2
p <- ggplot(pca_data, aes(x=Axis1, y=Axis2, colour=Pathotype,
                          # size=Pathogenicity,
                          shape=Sequencing_Platform)) + geom_point(size=4.5, alpha = 0.8) +
  # scale_color_gradient2() +
 scale_color_manual(values = path_cols) +
  # scale_size_manual(values = path_sizes[pca_data$Pathogenicity]) +
  geom_text_repel(aes(label=ifelse(Isolate %in% outliers, Isolate, ""))) +
  plot_theme()
# Plot and add labels (represent component variance as a percentage)
p +
  labs(shape="Seq. Platform", x=glue("C1: {round(percentVar[1]*100,2)}% variance"),
       y=glue("C1: {round(percentVar[2]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("{analysis_basename}_PCA_patho_platform"),
                ext = ".pdf", outdir = "plots",
                dateformat = FALSE), width = 10, height=7)


#### GenABEL ####

# Used SnpSift to convert from vcf to PLINK tped
files_prefix <- glue("./data/intermediate_files/{analysis_basename}.snpsift")
  convert.snp.tped(tpedfile = glue("{files_prefix}.tped"),
                   tfamfile = glue("{files_prefix}.tfam"),
                   out=glue("{files_prefix}.raw"),
                   strand="+")



pheno_data <- strata %>% rename(id=INDIVIDUALS) %>% mutate(sex=1) %>%
  .[c("id", "sex", "Pathogenicity")]
write_tsv(pheno_data, glue("{files_prefix}.dat"))
# Loads the snp and phenotypic data
snp_data <- load.gwaa.data(phenofile = glue("{files_prefix}.dat"),
                           genofile = glue("{files_prefix}.raw"),
                           force = TRUE, makemap = FALSE, sort = TRUE)
# compute IBS based on a random sample of 1000 autosomal marker
a <- ibs(snp_data)
a[1:5,1:5]
mds <- cmdscale(as.dist(1-a))
plot(mds)

qc1 <- check.marker(snp_data, maf = 0.01, fdrate = 0.01, ibs.threshold = 0.999, ibs.mrk = "all")
clean_snps <- snp_data[qc1$idok, qc1$snpok]
# Export clean geno table
geno_table <- as.data.frame(as.character(gtdata(clean_snps)))
# Fix Contig names
colnames(geno_table) <- gsub("Arab_Me14_", "", colnames(geno_table))
export_snps <- geno_table %>%
  tibble::rownames_to_column(var = "id") %>%
  inner_join(strata, ., by=c("INDIVIDUALS"="id")) %>% dplyr::rename(id=INDIVIDUALS, Patho_level=STRATA) %>%
  as.data.frame(.) %>% column_to_rownames("id")
row.names(export_snps) <- export_snps$id
# head(clean_snps@gtdata@snpnames)

#### SNPassoc analysis  ####
SNP.info.pos <- data.frame(snp=colnames(geno_table),
                           chr=sub("_\\d+", "", colnames(geno_table)),
                           pos=sub(".+_\\(d+)", "\\1", colnames(geno_table)))
myData <- setupSNP(data=export_snps,colSNPs=5:ncol(export_snps),sep="/", info=SNP.info.pos, sort = TRUE)
# Check missingness
# plotMissing(myData)
# Perform association test to Pathogenicity level
codom_res <- WGassociation(Pathogenicity, data=myData, model="codominant") # try different models
# Check summary
p_val <- 1e-3
min_snps <- 50
FDR <- 0.05

# Calculate qvalue
res_selected <- codom_res[!is.na(codom_res$codominant)]
# Select only chromosomes with more than 50 SNPs or ones with significantly associated ones

res_df <- cbind(attr(res_selected, "gen.info"), "pvalue"=res_selected$codominant) %>%
  bind_cols(., qvaluebh95(.$pvalue, FDR))
res_sum <- res_df %>% group_by(chr) %>% summarise(nSNPs=n(), significant=sum(significant))
chr_select <- res_sum %>% filter(nSNPs>=min_snps | significant>0) %>% droplevels() %>% .$chr
# Select SNPs in selected chromosomes
filt_snps <- res_df[res_df$chr %in% chr_select,]

signif_snps <- filt_snps %>% filter(significant==TRUE) %>% mutate(snp=as.character(snp)) %>% .$snp

# Save data
save.image(file = filedate("association_analysis_all", ".RData", "data"))
pdf(filedate("codom_SNPs", ".pdf", "plots"), width = 7, height = 6)
plot(res_selected, sort.chromosome=TRUE, centromere=FALSE, cutPval = c(0, p_val, 1), ylim.sup=1e-15)
dev.off()
# sign_snps <- names(codom_res$comments)[!is.na(codom_res$codominant) & codom_res$codominant<=p_val]
# Check association for a single SNP


# Perform stratified analysis (when we'll have more isolates in each region)
log_add_res <- scanWGassociation(Pathogenicity, data=myData, model="log-add")
