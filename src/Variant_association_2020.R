devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
# install.packages("devtools")
# install.packages("pacman")

# Install CRAN archived packages
devtools::install_version("GenABEL.data",version="1.0.0")
devtools::install_version("GenABEL",version="1.8-0")
devtools::install_bitbucket(repo = "bucklerlab/rtassel", ref = "master")
CRAN_packages <- c("tidyverse", "RColorBrewer", "ggrepel","GenABEL", "outliers", "adegenet", "SNPassoc", "vcfR", 
                   "glue", "paletteer", "here")
pacman::p_load(char=CRAN_packages)

pacman::p_load(tidyverse, here, scales, Microsoft365R, ISOweek,
               janitor, readxl)

# install.packages("SNPassoc")
# install.deps(CRAN_packages)

# install.deps("SNPassoc", repo = "bioc")
# Download working version of SNPassoc from http://www.creal.cat/media/upload/arxius/jr/SNPassoc/SNPassoc_1.8-5.zip and copy to R library

#### Read data files ####
# Read sample metadata and save to file
analysis_basename="A_rabiei_2018"
variant_method <- "snippy_multi"
analysis_folder <- "../Snippy_multi_17_07_2019/"
analysis_outdir <- glue("./output/{variant_method}/")
# read isolate information from excel file
haplotype_info <- readxl::read_excel("./sample_info/5_year_complete_SSR_db.xlsx", sheet = "Feb_2019") %>% 
  dplyr::select(Region, Isolate=Ind, State, Collection_Year=Year, Haplotype)
# read isolate sequencing info
sequencing_table <- readxl::read_excel("./sample_info/A_rabiei_isolate_list_for_wgs.xlsx",sheet = "submission_info") 
# create a dictionarty for sample-isolate 
sequencing_dict <- set_names(sequencing_table$Isolate, sequencing_table$Submission_id)

# load scoring table
scoring_table <- readxl::read_excel("./sample_info/A_rabiei_isolate_list_for_wgs.xlsx", "Path_score") %>% 
  filter(Host!="Seamer") %>% gather(key = "path_levels", value = "score", -Host) 
scoring_set <- unique(scoring_table$Host)
assign_score <- function(path_level, host, scoring_mat=scoring_table){
  scoring_mat %>% dplyr::filter(Host==host, path_levels %in% path_level) %>% dplyr::select(score) %>% as.numeric()
}

# define pathotypes
path_levels <- c("Low","Medium","Moderate","High","Very High", "Extreme") %>% factor(., levels = .)
  
path_groups <- setNames(paste0("Group", 0:5), levels(path_levels))
# Load sample and mating-type tables
mat_table <-  readxl::read_excel("./sample_info/2017_A_rabiei_DNA_for_MT_PCR.xlsx",sheet = "Sheet1") %>%
  dplyr::select(Isolate=Sample_name, MAT1_2=Amplified) %>%
  # left_join(readxl::read_excel("./sample_info/2017_A_rabiei_DNA_for_MT_PCR.xlsx",sheet = "11_04") %>%
  #             dplyr::select(Isolate=Sample_name, MAT1_2b=`MAT Amplified2`)) %>% 
  mutate(MAT_Type=if_else(MAT1_2=="Yes", "MAT1-2", "NA")) %>%
  dplyr::select(Isolate, MAT_Type)

samples_table <- readxl::read_excel("./sample_info/A_rabiei_isolate_list_for_wgs.xlsx", "Sequenced")  %>%
  dplyr::select(-Haplotype, Host=`Host Cultivar`, Year=Collection_Year) %>%
  mutate(Host=sub("PBA ", "", Host)) %>%
  left_join(mat_table) %>%
  left_join(haplotype_info %>% dplyr::select(Isolate, Haplotype)) %>%
  mutate_at(scoring_set, ~str_replace_na(., "N/A")) %>%
  mutate_if(is.character, ~str_replace_na(., "Unknown")) %>% mutate(Haplotype=str_replace(Haplotype, "Unknown", "TBD"))
# 
samples_table$Pathogenicity <- rowSums(map_dfc(scoring_set, function(n) map_dbl(samples_table[[n]], ~assign_score(., n))))
samples_table <- samples_table %>% dplyr::select(-one_of(scoring_set)) %>%
  mutate(Pathotype=path_levels[Pathogenicity+1], Patho.Group=paste0("Group",Pathogenicity)) %>%
  write_xlsx(.,
             "./sample_info/A_rabiei_isolate_list_for_wgs.xlsx",
             "sample_details_full", overwritesheet = TRUE)
samples_table <- readxl::read_excel("./sample_info/A_rabiei_isolate_list_for_wgs.xlsx",
                                    "sample_details_full")
  #,
# add read numbers from the server
read_stats <-   recent_file("raw_data/", ".+raw_reads.tsv") %>% read_tsv(., col_names = c("filename", "pair", "reads")) %>%
    mutate(Sample_id=sub("_R[12].f.+q", "", filename)) #%>% mutate(Isolate=sequencing_dict[Sample_id])

# Read mapping stats (after preparation by the markdown notebook)
# mapping_stats <- recent_file("raw_data/",
#                              glue("{analysis_basename}.+mapping.stats.txt")) %>% read_csv() %>%
#   inner_join(read_stats %>% group_by(Isolate) %>%
#                summarise(paired_reads=mean(reads))) %>%
#   write_xlsx(., glue("{analysis_outdir}/{analysis_basename}_mapping.xlsx"),
#                             "mapping_details")
mapping_stats_outfile <- recent_file("output/results//",
                            glue("A_rabiei_2018.+mapping.stats.xlsx"))
mapping_stats <- mapping_stats_outfile %>% 
  readxl::read_excel(., "sequencing_stats") %>% #  "mapping_details"
  inner_join(read_stats %>% group_by(Sample_id) %>%
               summarise(paired_reads=mean(reads)))

# summarise stats per batch
# mapping_summary <- recent_file(glue("{analysis_outdir}"), 
#                                glue("{analysis_basename}.+mapping.sum.txt")) %>%  read_tsv() %>% 
mapping_summary <- mapping_stats %>%  group_by(Sequencing_Centre) %>% 
              summarise(Coverage=sprintf("x%.2f", mean(Coverage)), 
                        Mapping_rate=mean(Mapping_rate), Mapping_qual=mean(Mapping_quality_mean), 
                        Total_paired_reads=sum(paired_reads)/1e6, Files=n()) %>% 
  mutate(Bases=if_else(Sequencing_Centre=="Macrogen",Total_paired_reads/10, 
                       Total_paired_reads*15/100)) %>% 
write_xlsx(., mapping_stats_outfile, 
           "batches_sum", overwritesheet = TRUE)

# combine stats from multiple sequencing samples/batches per isolate
sequenced_isolates <- mapping_stats %>% group_by(Isolate) %>% 
  summarise(GC_percentage=mean(GC_percentage), 
            Mapping_quality_mean=mean(Mapping_quality_mean), 
            Coverage=sum(Coverage)) %>% inner_join(samples_table, .) %>% 
  write_xlsx(.,  mapping_stats_outfile, 
             "mapping_per_sample")


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

# path_levels <- c("Low", "Medium", "Moderate", "High", "Very High")
# strata <- sequenced_isolates %>% # filter(!is.na(Sequencing_Centre)) %>%
#   mutate(INDIVIDUALS=Isolate, STRATA=path_levels[Pathogenicity+1],
#          Site=paste(Site, State, sep=", ")) %>%
#   dplyr::select(INDIVIDUALS, STRATA, Pathogenicity, Site, State, Collection_Year, Coverage) # %>%
#   # inner_join(mapping_per_isolate[c("Isolate", "Coverage")],
#   #                                     by= c("INDIVIDUALS"="Isolate"))  %>% 
#   write_tsv(filedate(analysis_basename, ".strata",glue("{analysis_outdir}/results"), dateformat = FALSE))



#### VCF EDA with vcfR ####
# vcf_file <- recent_file("./data", glue::glue("{analysis_basename}.+.vcf"))
vcf_file <- recent_file(analysis_folder, glue("{analysis_basename}_isolates_{variant_method}.Qual10.+poly.snps.vcf"))
vcf_file <- file.path(analysis_folder, "core.vcf")
vcf_basename <- tools::file_path_sans_ext(basename(vcf_file))
vcf <- read.vcfR(vcf_file)
# Which sequenced isolates are not in the vcf file
sequenced_isolates %>% filter(!Isolate %in% colnames(extract.gt(vcf))) %>% 
  dplyr::select(Isolate)

# Fix sample_names
# colnames(vcf@gt)[colnames(vcf@gt) %in% agrf_table$Submission_ID] <- agrf_dict[colnames(vcf@gt)[colnames(vcf@gt) %in% agrf_table$Submission_ID]]


#### Analysis ####
# filtered_vcf <- read.vcfR(filtered_vcf_file)
genind_obj <- vcfR2genind(vcf) # clean_vcf
# output as tfam format
sample_cols <- colnames(extract.gt(clean_vcf))
tfam_table <- samples_table %>% filter(Isolate %in% sample_cols) %>% 
  write_xlsx(., glue("{analysis_outdir}/results/{analysis_basename}_mapping.xlsx"), "vcf_samples", append=TRUE) %>% 
  dplyr::select(Host, Isolate, Pathogenicity) %>%
  mutate(V3=0, V4=0, V5=0, V6=Pathogenicity) %>%
  write_tsv(filedate(glue("{analysis_basename}_host_patho"), ".tfam", analysis_outdir, dateformat = FALSE), col_names = FALSE)

# load to adegenet genlight (straight from tped file)
# obj <- genomic_converter("data/A_rabiei_on_me14_DP_GQ_corr.bt2.fb.vcf", output = "genind", strata = paste0(analysis_basename, ".strata"),  snp.ld = "maf",parallel.core = 0) #  imputation.method="rf",
# # fix chromosome names
# obj$tidier_data <- obj$tidy.data %>% mutate(MARKERS=sub("Arab_Me14_", "", MARKERS), CHROM=sub("Arab_Me14_", "", CHROM))

# Analyze missingness
# ibm.rabiei <- missing_visualization("data/A_rabiei_on_me14.bt2.fb.vcf", strata = "data/A_rabiei_isolates.strata", strata.select = c("STRATA", "Site"))
# genind_obj <- obj$genind.no.imputation

#### Perform some EDA  ####
samples_factors <- samples_table %>% mutate(State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA")),
  # Year=factor(Year, levels=as.character(min(Year, na.rm = TRUE):max(Year, na.rm = TRUE))),
  # Host=fct_relevel(fct_infreq(fct_lump_min(Host, min=6)), "Other", after = Inf),
  # Haplotype=fct_relevel(fct_infreq(fct_lump_min(Haplotype, min=2)), "Other", after = Inf),
  Patho.Group=factor(Patho.Group, levels=paste0("Group", 0:5))) %>%
  mutate_at(vars(Host, Haplotype, Patho.Group), ~fct_explicit_na(., na_level = "Unknown"))

# Convert to allele frequencies (manage missing data)
X <- adegenet::tab(genind_obj, freq = TRUE, NA.method = "mean")
pca1 <- ade4::dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 20)

pca_data <- pca1$li %>% rownames_to_column("Isolate")  %>%
  inner_join(samples_factors) %>% # sequenced_isolates
  mutate(seq_fontface=if_else(Sequenced=="Illumina", "plain", "bold")) %>% 
         # Pathotype=factor(Pathotype, levels=path_levels)) %>% #,
         # Coverage=round(Coverage, 2)) %>%
  left_join(sequencing_table %>% filter(Sequencing_Centre=="Macrogen") %>% 
  dplyr::select(Isolate, Sequencing_Centre)) %>% 
  mutate(Sequencing_Centre=if_else(is.na(Sequencing_Centre), "AG", Sequencing_Centre)) %>% 
  left_join(sequenced_isolates %>% dplyr::select(Isolate, GC_percentage, Mapping_quality_mean, Coverage )) %>% 
  arrange(Year)

  # mutate(Pathotype=factor(STRATA, levels=path_levels)Host)
  #        
  #        ,
  #     Sequencing_Platform=if_else(Sequencing_Centre=="AGRF", "NextSeq", "HiSeq2500")) %>%
  # filter(!is.na(Pathotype))

# Calculate variance for each component (columns starting with 'C')
pca_var <- pca1$li %>% summarise_at(vars(starts_with("Axis")), var)
# Calculate the percentage of each component of the total variance
percentVar <- pca_var/sum(pca_var)
# Define colour palette, but get rid of the awful yellow - number 6
my_colours = list(
  State = levels(pca_data$State) %>% setNames(as.character(paletteer_d("RColorBrewer::Set1", length(.))), .),
  Year = sort(unique(pca_data$Year)) %>% setNames(as.character(paletteer_d("ggsci::default_uchicago", length(.))), .),
  Patho.Group = levels(pca_data$Patho.Group) %>% 
    setNames(as.character(paletteer_d("RColorBrewer::RdYlGn", direction = -1))[round(seq(4,11, length.out = length(.)) ,0)], .),
  Haplotype = levels(pca_data$Haplotype) %>% setNames(as.character(paletteer_d("RColorBrewer::Dark2", length(.))), .)
)


pal <- brewer.pal(9, "Set1")
pal2 <- brewer.pal(8, "Dark2")


# path_guide <- tibble(path_levels, col=pal[c(3,2,5,1,4)], size=seq(2,4,length.out = 5))
# path_sizes <- setNames(seq(2,4,length.out = 5), path_levels)
path_cols <- my_colours$Patho.Group# setNames(pal[c(3,2,5,1,4,9)], path_levels)
seq_cols <- adjustcolor( c("grey15", "dodgerblue3"), alpha.f = 0.8)
seq_face <- setNames(c("bold", "plain"),unique(pca_data$Sequenced))

shapes <- c(21:25)
# Isolate to annotate

# outliers <- pca_data %>% filter(grepl("PacBio", Sequenced, ignore.case = TRUE)) %>% .[,"Isolate"]
# outliers <- c("16CUR017", "TR9543", "TR9529", "15CUR003") # PacBio-sequenced
# outliers <- pca_data$Isolate
#### PCA plot 1-2 biological ####
# shape - Year, Seq.Provider, Host (State)
# Colors - Haplotype (diverged), State (diverged), Pathotype (RdYlGn)
# Size - Pathogenicity

# Create the plot for C1 and C2
# Create the plot for C1 and C2
plot_comps <- 1:2
plot_axes <- paste0("Axis", plot_comps)
outliers <- pca_data %>% filter(abs(!!sym(plot_axes[2])) > 2.5 | abs(!!sym(plot_axes[1])) > 2.5) %>% .[,"Isolate"]
ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), fill=Patho.Group, shape=Host)) + 
# p <- ggplot(pca_data, aes(x=Axis1, y=Axis2, fill=Pathotype, 
#                           # size=Year,
#                           shape=Host)) # Host
# Plot and add labels (represent component variance as a percentage)
  geom_point(alpha = 0.8, size=5) + # , stroke = 1
  # geom_text_repel(aes(colour=Sequenced, label=Isolate, fontface=seq_fontface),
  # show.legend = FALSE,size=3.5, point.padding = 0.75) +
  geom_text_repel(aes(label=if_else(Isolate %in% outliers, glue("{Isolate} ({Year})"), ""),
                      colour=Sequenced, fontface=seq_fontface),
                      # colour=adjustcolor( c("dodgerblue3"), alpha.f = 0.8)),
              show.legend = FALSE, size=3.5, point.padding = 0.75) +
  # scale_color_() +
  scale_fill_manual(values = my_colours$Patho.Group) + 
  # scale_fill_brewer(palette = "RdYlGn", direction = -1) + 
  # scale_fill_manual(values = pal, aesthetics = c("fill")) + 
  scale_shape_manual(values = shapes) +
  scale_color_manual(values=adjustcolor( seq_cols, alpha.f = 0.8)) + # [1]
  # scale_size_continuous(range = c(4,7)) +
  guides(shape = guide_legend(override.aes = list(size = 5, fill=pal[2]), order = 1),
         fill = guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2) ,
         size = guide_legend(override.aes = list(fill=pal[2], pch=shapes[1]), order = 3)) +
  #  +
  # scale_size_manual(values = path_sizes[pca_data$Pathogenicity]) +
  
  # geom_text_repel(aes(label=Haplotype), size=4, point.padding = 0.75) +
  plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue("C{plot_comps[1]}: {round(percentVar[plot_comps[1]]*100,2)}% variance"),
       y=glue("C{plot_comps[2]}: {round(percentVar[plot_comps[2]]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("{analysis_basename}_PCA_hostShape_pathoFill_labelOutliers_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = glue("{analysis_outdir}/plots")), width = 8.5, height=6.5)
#### PCA plot 1-2 biological2 ####
# shape - Year, Seq.Provider, Host (State)
# Colors - Haplotype (diverged), State (diverged), Pathotype (RdYlGn)
# Size - Pathogenicity
# Create the plot for C1 and C2
plot_comps <- 3:4
plot_axes <- paste0("Axis", plot_comps)
p <- ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), fill=Haplotype, 
                          size=Year,
                          shape=State)) # Host
# Plot and add labels (represent component variance as a percentage)
p + geom_point(alpha = 0.8) + # , size=5, stroke = 1
  # geom_text_repel(aes(colour=Sequenced, label=Isolate, fontface=seq_fontface), show.legend = FALSE,size=3.5, point.padding = 0.75) +
  # geom_text_repel(aes(label=if_else(Isolate %in% outliers, glue("{Isolate} ({Year})"), ""), colour=adjustcolor( c("dodgerblue3"), alpha.f = 0.8)),
  #                 show.legend = FALSE, size=3.5, point.padding = 0.75) +
  # scale_color_() +
  # scale_fill_brewer(palette = "Dark2", direction = 1) + 
  scale_fill_manual(values = my_colours$Haplotype) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values=adjustcolor( seq_cols[1], alpha.f = 0.8)) +
  scale_size_continuous(range = c(3.5,6.5)) +
  guides(shape = guide_legend(override.aes = list(size = 5, fill=pal[2]), order = 1),
         fill = guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2) ,
         size = guide_legend(override.aes = list(fill=pal[2], pch=shapes[1]), order = 3)) +
  #  +
  # scale_size_manual(values = path_sizes[pca_data$Pathogenicity]) +
  
  # geom_text_repel(aes(label=Haplotype), size=4, point.padding = 0.75) +
  plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue("C{plot_comps[1]}: {round(percentVar[plot_comps[1]]*100,2)}% variance"),
        y=glue("C{plot_comps[2]}: {round(percentVar[plot_comps[2]]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("{analysis_basename}_PCA_yearSize_stateShape_haploFill_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = glue("{analysis_outdir}/plots")), width = 9, height=7)
#### PC3-4 plot ####
# Isolate to annotate
outliers <- pca_data %>% filter(abs(Axis4) > 1 | abs(Axis3)>1) %>% .[,"Isolate"]
# Create the plot for C1 and C2
p <- ggplot(pca_data, aes(x=Axis3, y=Axis4, fill= State, 
                          # size=Year,
                          shape=Host)) # Host
# Plot and add labels (represent component variance as a percentage)
p + geom_point(alpha = 0.8, size=5) + # , stroke = 1
  # geom_text_repel(aes(colour=Sequenced, label=if_else(Isolate %in% outliers, Isolate, ""), fontface=seq_fontface), show.legend = FALSE,size=3.5, point.padding = 0.75) +
  geom_text_repel(aes(label=if_else(Isolate %in% outliers, glue("{Isolate} ({Year})"), ""), 
                      colour=adjustcolor( c("dodgerblue3"), alpha.f = 0.8)),
              show.legend = FALSE, size=3.5, point.padding = 0.75) +
  # scale_color_() +
  scale_fill_brewer(palette = "Set1", direction = 1) + 
  # scale_fill_manual(values = pal, aesthetics = c("fill")) + 
  scale_shape_manual(values = shapes) +
  scale_color_manual(values=adjustcolor(seq_cols[1], alpha.f = 0.8)) +
  # scale_size_continuous(range = c(3,6)) +
  guides(shape = guide_legend(override.aes = list(size = 5, fill=path_cols[2]), order = 1),
         fill = guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2) ,size = guide_legend(override.aes = list(fill=path_cols[2], pch=shapes[1]), order = 3)) +
  #  +
  # scale_size_manual(values = path_sizes[pca_data$Pathogenicity]) +
  
  # geom_text_repel(aes(label=Haplotype), size=4, point.padding = 0.75) +
  plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue("C3: {round(percentVar[3]*100,2)}% variance"),
       y=glue("C4: {round(percentVar[4]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("{analysis_basename}_PCA_host_state_labelsyear_PC3-4"),
                ext = ".pdf", outdir = glue("{analysis_outdir}/plots"),
                dateformat = FALSE), width = 8.5, height=7)


#### PCA plot 1-2 technical ####
# shape - Year, Seq.Provider, Host (State)
# Colors - Haplotype (diverged), State (diverged), Pathotype (RdYlGn)
# Size - Pathogenicity
# Create the plot for C1 and C2
plot_comps <- 3:4
plot_axes <- paste0("Axis", plot_comps)
p <- ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]),
                          fill=GC_percentage,
                          size=Coverage,
                          shape=Sequencing_Centre))#, fill=Sequencing_Centre)) # Host
# Plot and add labels (represent component variance as a percentage)
p + geom_point(alpha = 0.8)+ #, size=5) + # , stroke = 1
  # geom_text_repel(aes(colour=Sequenced, label=Isolate, fontface=seq_fontface), show.legend = FALSE,size=3.5, point.padding = 0.75) +
  # geom_text_repel(aes(label=Coverage, colour=adjustcolor( c("dodgerblue3"), alpha.f = 0.8)), 
  #             show.legend = FALSE, size=3.5, point.padding = 0.75) +
  # scale_color_() +
  # scale_fill_brewer(palette = "Blues", direction = 1) + 
  scale_fill_paletteer_c("viridis::inferno")+
  # scale_fill_manual(values = pal, aesthetics = c("fill")) + 
  scale_shape_manual(values = shapes) +
  scale_color_manual(values=seq_cols) +
  # scale_size_continuous(range = c(3.5,6.5)) +
  guides(shape = guide_legend(override.aes = list(size = 5, fill=pal[2]), order = 1, title = "Seq.Centre"),
         fill = guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2, title = "GC (%)") ,
         size = guide_legend(override.aes = list(fill=pal[2], pch=shapes[1]), order = 3)) +
  #  +
  # scale_size_manual(values = path_sizes[pca_data$Pathogenicity]) +
  
  # geom_text_repel(aes(label=Haplotype), size=4, point.padding = 0.75) +
  plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue("C{plot_comps[1]}: {round(percentVar[plot_comps[1]]*100,2)}% variance"),
        y=glue("C{plot_comps[2]}: {round(percentVar[plot_comps[2]]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("{analysis_basename}_PCA_technical_viridis_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = glue("{analysis_outdir}/plots")), width = 8.5, height=7)

save.image(filedate(glue("{analysis_basename}_PCA_plots"), ext = ".RData", glue("{analysis_outdir}")))
#### GenABEL ####

# Used SnpSift to convert from vcf to PLINK tped
# files_prefix <- glue("./data/intermediate_files/{analysis_basename}.snpsift")
#   convert.snp.tped(tpedfile = glue("{files_prefix}.tped"),
#                    tfamfile = glue("{files_prefix}.tfam"),
#                    out=glue("{files_prefix}.raw"),
#                    strand="+")
# 
# 
# 
# pheno_data <- strata %>% rename(id=INDIVIDUALS) %>% mutate(sex=1) %>%
#   .[c("id", "sex", "Pathogenicity")]
# write_tsv(pheno_data, glue("{files_prefix}.dat"))
# # Loads the snp and phenotypic data
# snp_data <- load.gwaa.data(phenofile = glue("{files_prefix}.dat"),
#                            genofile = glue("{files_prefix}.raw"),
#                            force = TRUE, makemap = FALSE, sort = TRUE)
# # compute IBS based on a random sample of 1000 autosomal marker
# a <- ibs(snp_data)
# a[1:5,1:5]
# mds <- cmdscale(as.dist(1-a))
# plot(mds)
# 
# qc1 <- check.marker(snp_data, maf = 0.01, fdrate = 0.01, ibs.threshold = 0.999, ibs.mrk = "all")
# clean_snps <- snp_data[qc1$idok, qc1$snpok]

# head(clean_snps@gtdata@snpnames)

#### GAPIT analysis ####
BiocManager::install("snpStats")
pacman::p_load(char = c("multtest", "gplots", "LDheatmap", "genetics", "EMMREML", 
                        "scatterplot3d", "compiler"))
pacman::p_load_gh("jiabowang/GAPIT3")
# source("http://zzlab.net/GAPIT/GAPIT.library.R")
# source("http://zzlab.net/GAPIT/gapit_functions.txt")
# source("http://zzlab.net/GAPIT/emma.txt")
options(stringsAsFactors = FALSE)

gapit_phenotypes <- read_tsv("sample_info/230_GAPIT_phenotypes.txt")


# prepare phenotypic data
pheno_data <- samples_table %>% as.data.frame(., stringsAsFactors = FALSE) %>% 
  dplyr::select(Taxa=Isolate, Pathogenicity, Year) # , Host, Haplotype, State

# Create HapMap genotypic file
clean_vcf <- vcf
GT_data <- extract.gt(clean_vcf, return.alleles = TRUE) %>% data.frame(.) %>% 
  setNames(., colnames(clean_vcf@gt)[-1])
# GT_data[GT_data==1] <- 2
GT_data[is.na(GT_data)] <- "N"
GT_data[1:10,1:15]  
# SNP info
dummy_data <-  replicate(7, rbind(rep(NA_character_, nrow(GT_data))), simplify = "matrix") %>% 
  as.data.frame(.)
snp_info <-   getFIX(clean_vcf) %>% as.data.frame(.) %>% 
  dplyr::mutate(chrom=CHROM, pos=POS, rs=paste(CHROM, POS, sep = "_"), 
                                           alleles=paste(REF, ALT, sep = "/") ) %>% 
  dplyr::select(rs,  alleles, chrom, pos) 
GT_hapmap <- dplyr::bind_cols(snp_info, dummy_data, GT_data) %>% rbind(colnames(.), .)
 
# Create numeric genotypic file
# GD_data <- t(extract.gt(clean_vcf, as.numeric = TRUE)) %>% as.data.frame(.)
# GD_data[GD_data==1] <- 2
# # GD_data[is.na(GD_data)] <- NA
# 
# # SNP info
# GM_data <- data.frame(Name=colnames(GD_data), 
#                         Chromosome=getCHROM(clean_vcf), Position=getPOS(clean_vcf)) 
# 
# GD_data <- GD_data  %>% tibble::rownames_to_column(var = "taxa") 


#Step 2: Run GAPIT
myGAPIT <- GAPIT(Y=pheno_data,G = GT_hapmap ,PCA.total=8, Model.selection = TRUE, 
                 model = c("CMLM", "GLM","MLM","MLMM","FarmCPU", "SUPER"),
                 NJtree.group=4,                                       # set the number of clusting group in Njtree plot
                 # QTN.position=mysimulation$QTN.position,
                 Inter.Plot=TRUE,                                      # perform interactive plot
                 Multiple_analysis=TRUE,                               # perform multiple analysis
                 PCA.3d=TRUE,                                          # plot 3d interactive PCA
                 file.output=T, )





  
#### SNPassoc analysis  ####
# Export clean geno table
geno_table <- extract.gt(clean_vcf, as.numeric = TRUE)
geno_table[geno_table==1] <- 2
snp_table <- apply(geno_table, 1, function(s) snp(s,name.genotypes=0:2))
snp_table[1:10,1:10]

# Fix Contig names
colnames(snp_table) <- gsub("Arab_Me14_", "", colnames(snp_table))
export_snps <- snp_table %>% as.data.frame(.) %>% 
  tibble::rownames_to_column(var = "Isolate") 
  inner_join(sequenced_isolates, .) %>% # dplyr::rename(id=INDIVIDUALS, Patho_level=STRATA) %>%
   column_to_rownames("Isolate")
export_snps[1:10,1:10]  
# row.names(export_snps) <- export_snps$id
# example <- data(SNPs)
# example_info <- data(SNPs.info.pos)
SNP.info.pos <- data.frame(snp=colnames(snp_table),
                           chr=sub("_\\d+", "", colnames(snp_table)),
                           pos=sub(".+_(\\d+)", "\\1", colnames(snp_table)))
SNP.info.pos[1:10,]
SNP.info.pos %>% group_by(chr) %>% summarise(n()) %>% print(n=Inf)
myData <- setupSNP(data=export_snps,colSNPs=ncol(sequenced_isolates):ncol(export_snps),sep="/", 
                   info=SNP.info.pos, sort = TRUE)
# Check missingness
# plotMissing(myData)
# Perform association test to Pathogenicity level
res <- WGassociation(Pathogenicity, data=myData, model="all") # try different models
# res <- scanWGassociation(Pathogenicity, data=myData, model="all")
# Check summary
p_val <- 1e-3
min_snps <- 50
FDR <- 0.05

# Calculate qvalue
res_selected <-res[!is.na(res$`log-additive`)]
# Select only chromosomes with more than 50 SNPs or ones with significantly associated ones

res_df <- cbind(attr(res_selected, "gen.info"), "pvalue"=res_selected$`log-additive`) %>%
  bind_cols(., qvaluebh95(.$pvalue, FDR))
# check if it found ny significant SNPs
res_df %>% filter(significant!=FALSE)
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
