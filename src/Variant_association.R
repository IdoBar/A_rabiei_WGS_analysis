devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# Install these to enable installation of radiator and grur
# install_cran(c("future", "Rcpp", "apex", "copula", "gsl", "ADGofTest", "stabledist", "pcaPP", "pspline", "shinyFiles",
#                "RJSONIO", "swfscMisc", "mapdata", "maps"))
# install.deps(c("tidyverse/glue", "tidyverse/tidyselect", "thierrygosselin/radiator", "thierrygosselin/grur"), repo="git") # need to be run first, otherwise fails to re-install tidyverse packages if they're already loaded
install.deps(c("thierrygosselin/radiator"), repo="git")
CRAN_packages <- c("tidyverse", "RColorBrewer", "ggrepel","GenABEL", "adegenet", "SNPassoc", "here")
# install.packages("SNPassoc")
pacman::p_load(char = CRAN_packages)

# install.deps("SNPassoc", repo = "bioc")
# Download working version of SNPassoc from http://www.creal.cat/media/upload/arxius/jr/SNPassoc/SNPassoc_1.8-5.zip and copy to R library
# # Theme setup
# plot_theme <-  theme_grey(base_size=20) +
#   theme(axis.title.y=element_text(face="bold", vjust = 1.5, size=rel(0.8)),
#         axis.title.x=element_text(face="bold", vjust = 0.1, size=rel(0.8)),
#         legend.title=element_text(size=rel(0.8), hjust=0.5),
#         legend.text=element_text(size = rel(0.7),lineheight = 1.5),
#         #panel.grid.minor=element_blank(),
#         strip.text=element_text(size=rel(0.6)),
#         strip.switch.pad.grid = unit(0, "lines"),
#         strip.background=element_rect(fill = "lightskyblue"))#, colour = "black", size

# Read sample metadata and save to file
samples_table <- readxl::read_excel("data/P_rabiei_isolate_list_sent_for_sequencing.xlsx")
tfam_table <- samples_table %>% filter(Sequenced=="Yes") %>% dplyr::select(State, Isolate) %>%
  mutate(V3=0, V4=0, V5=0, V6=-9)
path_levels <- c("Low", "Medium", "Moderate", "High", "Very High")
strata <- samples_table %>% filter(Sequenced=="Yes") %>%
  mutate(INDIVIDUALS=Isolate, STRATA=path_levels[Pathogenicity+1], Site=paste(Site, State, sep="_")) %>%
  dplyr::select(INDIVIDUALS, STRATA, Pathogenicity, Site, State)
analysis_basename <- "data/A_rabiei_isolates"
if (!file.exists(paste0(analysis_basename, ".strata"))) write_tsv(strata, "data/A_rabiei_isolates.strata")


# load to adegenet genlight (straight from tped file)
obj <- genomic_converter("data/A_rabiei_on_me14.bt2.fb.vcf", output = "genind", strata = paste0(analysis_basename, ".strata"), imputation.method="rf", snp.ld = "random")
# fix chromosome names
obj$tidier_data <- obj$tidy.data %>% mutate(MARKERS=sub("Arab_Me14_", "", MARKERS), CHROM=sub("Arab_Me14_", "", CHROM))

# Analyze missingness
# ibm.rabiei <- missing_visualization("data/A_rabiei_on_me14.bt2.fb.vcf", strata = "data/A_rabiei_isolates.strata", strata.select = c("STRATA", "Site"))
# genind_obj <- obj$genind.no.imputation

#### Perform some EDA  ####
# Convert to allele frequencies (manage missing data)
X <- tab(genind_obj, freq = TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = TRUE)
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 4)

pca_data <- pca1$li %>% rownames_to_column(.)  %>% left_join(., strata, by = c( "rowname"= "INDIVIDUALS")) %>%
  mutate(Pathogenocity=factor(STRATA, levels=path_levels))

# Calculate variance for each component (columns starting with 'C')
pca_var <- pca1$li %>% summarise_at(vars(starts_with("Axis")), var)
# Calculate the percentage of each component of the total variance
percentVar <- pca_var/sum(pca_var)
# Define colour palette, but get rid of the awful yellow - number 6
pal <- brewer.pal(9, "Set1")[-6]
# Create the plot for C1 and C2
p <- ggplot(pca_data, aes(x=Axis1, y=Axis2, colour=Site, size=Pathogenicity)) + geom_point(size=3) +
 scale_color_manual(values = pal) +
  geom_text_repel(aes(label=rowname)) + plot_theme
# Plot and add labels (represent component variance as a percentage)
p + xlab(paste0("C1: ",round(percentVar[1]*100,2),"% variance")) +
  ylab(paste0("C2: ",round(percentVar[2]*100, 2),"% variance"))
# Save plot to a pdf file
ggsave(filedate("A_rabiei_isolates_PCA", ext = ".pdf", outdir = "plots"), width = 8, height=6)


#### GenABEL ####

# Used SnpSift to convert from vcf to PLINK tped
if (!file.exists(paste0(analysis_basename, ".raw"))) {
  convert.snp.tped(tpedfile = paste0(analysis_basename, ".tped"),
                   tfamfile = paste0(analysis_basename, ".tfam"),
                   out=paste0(analysis_basename, ".raw"),
                   strand="+")
}

pheno_data <- strata %>% rename(id=INDIVIDUALS) %>% mutate(sex=1) %>%
  .[c("id", "sex", "Pathogenicity")]
if (!file.exists(paste0(analysis_basename, ".dat"))) write_tsv(pheno_data, "data/A_rabiei_isolates.dat")
# Loads the snp and phenotypic data
snp_data <- load.gwaa.data(phenofile = paste0(analysis_basename, ".dat"),
                           genofile = paste0(analysis_basename, ".raw"),
                           force = TRUE, makemap = FALSE, sort = TRUE)
qc1 <- check.marker(snp_data, maf = 0.05, fdrate = 0.05)
clean_snps <- snp_data[qc1$idok, qc1$snpok]
export_snps <- as.data.frame(as.character(gtdata(clean_snps))) %>% tibble::rownames_to_column(var = "id") %>%
  left_join(strata, ., by=c("INDIVIDUALS"="id")) %>% rename(id=INDIVIDUALS, Patho_level=STRATA) %>%
  as.data.frame(.)
row.names(export_snps) <- export_snps$id
# head(clean_snps@gtdata@snpnames)
#### SNPassoc analysis  ####
SNP.info.pos <- data.frame(snp=clean_snps@gtdata@snpnames,
                           chr=sub("_\\d+", "", clean_snps@gtdata@snpnames), pos=clean_snps@gtdata@map)
myData<-setupSNP(data=export_snps,colSNPs=6:ncol(export_snps),sep="/", info=SNP.info.pos, sort = TRUE)
# Check missingness
# plotMissing(myData)
# Perform association test to Pathogenicity level
codom_res <- scanWGassociation(Pathogenicity, data=myData, model="codominant") # try different models
# Check summary
p_val <- 1e-5
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
save.image(file = "data/association_analysis.RData")
pdf(filedate("codom_SNPs", ".pdf", "plots"), width = 7, height = 6)
plot(res_selected, sort.chromosome=TRUE, centromere=FALSE, cutPval = c(0, p_val, 1), ylim.sup=1e-15)
dev.off()
# sign_snps <- names(codom_res$comments)[!is.na(codom_res$codominant) & codom_res$codominant<=p_val]
# Check association for a single SNP


# Perform stratified analysis (when we'll have more isolates in each region)
log_add_res <- scanWGassociation(Pathogenicity, data=myData, model="log-add")
