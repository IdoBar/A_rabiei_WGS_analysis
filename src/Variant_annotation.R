devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# devtools::install_github("jokergoo/ComplexHeatmap")
CRAN_packages <- c("tidyverse", "circlize", "RColorBrewer", "gplots",  "seqinr", "gridExtra", "grid",
                   "vcfR", "scales","ggrepel") # , "reutils" , "ComplexHeatmap", "GenABEL",
# pacman::p_load(char = CRAN_packages, update = TRUE)

pacman::p_load(char = CRAN_packages) # , force=TRUE

bioc_packages <- c("GenomicFeatures", "VariantAnnotation", "GenomicRanges", "Rsamtools")
                    #, "genphen", "GENESIS", "GWASTools")
# pacman::p_load(char = bioc_packages)
pacman::p_load(char = bioc_packages)

options(stringsAsFactors = FALSE)


# https://www.biostars.org/p/18954/
# Compress and index vcf file
vcf_file <- "../Snippy_multi_12_08_2019/core.vcf"
gff_file <- "../../Reference_genomes/Arab_me14_Augustus_genes.gff3"
genome_fa <- "../../Reference_genomes/Arab_me14.fasta"
vcf_head <- read_tsv(vcf_file, comment = "##", n_max = 10)
vcf_samples <- colnames(vcf_head)[10:ncol(vcf_head)]


#### Load and filter with VariantAnnotation ####
idx_vcf <- paste(vcf_file, "bgz", "tbi", sep = ".")
# zipped_vcf <- bgzip(vcf_file, paste0(vcf_file, ".bgz"))
idx <- indexTabix(bgzip(vcf_file, paste0(vcf_file, ".bgz"), overwrite = TRUE), "vcf")
core_vcf <- readVcf(idx)

# Filter externally with SnpSift and upload here, much more efficient
# Create filters
# call_rate <- 0.9
# min_called_samples <- floor(call_rate*length(vcf_samples))
# # Prefilters
# isSNP <- function(x) {
#   grepl("TYPE=snp", x, fixed=TRUE)
# }
# # Filters
# SampNum <- function(x, numPass=min_called_samples){
#   # x <- snp_vcf[39,]
#   gt <- geno(x)$GT
#   ns <- rowSums(!(is.na(gt) | gt=="." | gt=="./."))
#   ns>=numPass
# }
# 
# ReadDepth <- function(x, minDP=5, maxDP=1000, numPass=min_called_samples){
#   dp <- geno(x)$DP
#   test <- dp>minDP & dp<maxDP
#   total_dp <- rowSums(!is.na(test) & test)
#   total_dp >= numPass
# }
# # Check if marker is polymorphic
# PolyMorph <- function(x, hom_rate=0.05, biallelic=TRUE){
#   gt <- geno(x)$GT
#   poly <- apply(gt, 1, function(g) {
#     # props <- prop.table(table(g))
#     prop_names <- unique(g) # names(props)
#     log_vec <- ifelse(biallelic, length(prop_names[!grepl("\\.", prop_names)])==2,
#                       length(prop_names[!grepl("\\.", prop_names)])>1)
#     # p_table <- prop.table(table(as.character(g)))
#     # if (sum(grepl("[AB]", names(p_table)))==2) {
#     #   res <- p_table["A"]>hom_rate && p_table["B"]>hom_rate
#     # }  else res <- FALSE
#     # return(res)
#     
#   } )
# 
# }
# 
# # Add filters to rules
# filters <- FilterRules(list(SampNum, isSNV, PolyMorph))
# # Create a temp file for the filtered vcf
# # destination.file <- tempfile()
# destination.file <- sub(".vcf", ".PolyBi.snps.vcf", vcf_file, fixed = TRUE)
# filterVcf(paste0(vcf_file, ".bgz"),  "Arab_me14", destination.file, prefilters=FilterRules(isSNP),
#            param = ScanVcfParam(samples = vcf_samples),filters = filters, 
#           verbose=TRUE ) #
# snp_vcf <- readVcf(destination.file)
# samples(scanVcfHeader(destination.file))
# Or if filtering was performed before
snp_vcf <- core_vcf

# Change names to match fasta and txdb (remove prefix and leave only ctg##)
# seqlevels(snp_vcf) <- sub("Arab_Me14_", "", seqlevels(snp_vcf))

#rd <- rowRanges(snp_vcf) # extract the vcf content as a Grange



#### Load gene models  ####
txdb <- makeTxDbFromGFF(gff_file, format="gff3", dataSource = "CCDM Curtin",
                        organism = "Ascochyta rabiei",
                        taxonomyId = 5454)
# seqlevels(txdb) <- sub("Arab_Me14_", "", seqlevels(txdb))
# Break the vcf by each isolate
# Get a dataframe with all types of variants
getAllVariants <- function(vcf, txdb){
  coding_var <- locateVariants(vcf, txdb, CodingVariants())
  intron_var <- locateVariants(vcf, txdb, IntronVariants())
  idx <- GenomicRanges:::get_out_of_bound_index(vcf)
  intergenic_var <- locateVariants(vcf, txdb, IntergenicVariants())
  promoter_var <- locateVariants(vcf, txdb, PromoterVariants())
  allvars <- c(coding_var, intron_var, intergenic_var,promoter_var )  %>% trim(.) 
  return(as.data.frame(allvars,row.names = 1:length(allvars)))
}
# Find all variants for each isolate and add it up to a list of data frames
sample_names <- samples(header(snp_vcf))
varsDF <- tibble()
for (s in sample_names[1]){
  varsDF <- getAllVariants(snp_vcf[,s], txdb) %>% as_tibble() %>%  bind_rows(varsDF) %>% 
    distinct(seqnames,start,end, strand, LOCATION, LOCSTART,                                                                                                    LOCEND, QUERYID,  GENEID, .keep_all = TRUE) %>% mutate(CDSID=map_chr(CDSID, ~paste0(unlist(.x), collapse = ";")))
  LogMsg(glue::glue("Processed sample {s}, varsDF currently contains {nrow(varsDF)} unique variants and its memory allocation is estimated at {format(object.size(varsDF), units = 'auto')}"))
}
# varsDF %>% filter(GENEID=="gene3217")
# format(object.size(varsDF), units = 'auto')
# vars_by_sample <- sample_names[1:5] %>% purrr::map(~getAllVariants(snp_vcf[,.x], txdb))
# vars_DFlist <- DataFrameList(vars_by_sample)
# rm(vars_DFlist, varsDF)
# 
# # Combine tables
# varsDF <- as.data.frame(stack(vars_DFlist))
# Plot SNP density
x_ticks <- scales::cbreaks(c(0,3e6), breaks = seq(0,3e6, by = 1e6))
snpDensity <- ggplot(varsDF) +
  geom_histogram(aes(x=start, fill=LOCATION),binwidth=1e5) + # pick a binwidth that is not too small
  facet_wrap(~ seqnames,ncol=3)  + plot_theme("grey", 15) +
  scale_x_continuous(breaks=x_ticks$breaks, labels = x_ticks$labels) # + xlim(c(0,2e6)) + ylim(c(0,400))

ggsave(filedate("SNPDensity_medium", ".pdf", "output/snippy_multi/plots"), snpDensity,width = 10, height = 15)

# Load genome file to predict coding regions
# genome_fa <- "../A_rabiei_me14_short_names.fasta"
if (!file.exists(paste0(genome_fa, ".fai"))) indexFa(genome_fa)
genome_seqs <- read.fasta(genome_fa)
coding <- predictCoding(snp_vcf, txdb, FaFile(genome_fa)) %>% trim(.) 

codingDF <- as.data.frame(coding, row.names = 1:length(coding)) %>% mutate(varname=sub("(ctg\\d+):(\\d+)_.+", "\\1.\\2", names(coding)))

# combine association and annotaion results
# import significant SNPs from DArT analysis
pvalue_thresh <- 0.05
snp_flanking_region <- 5000
signif_snps <- readxl::read_excel("../../A_rabiei_DArT/output/A_rabiei_pathogenicity_snpassoc_results.xlsx") %>% filter(pvalue<=pvalue_thresh)
# find SNPs within the specified regions from the associated snps
assoc_snp_annot <- tibble()
for (i in 1:nrow(signif_snps)){
  assoc_snp_annot <- varsDF %>% mutate(chr=sub("Arab_Me14_", "", seqnames), varname=paste(seqnames, start, sep="."), gene_id=as.character(GENEID)) %>% 
    filter(chr %in% sub("^RYYQ0.+_(ctg\\d+)", "\\1",signif_snps$chr[i]),  start>=signif_snps$pos[i]-snp_flanking_region, # LOCATION!="intergenic",
           start<=signif_snps$pos[i]+snp_flanking_region) %>%   
    left_join(., codingDF[c("varname", "CONSEQUENCE", "REFAA", "VARAA")]) %>% 
    # dplyr::select(-name) %>% 
    bind_rows(assoc_snp_annot) %>%
    distinct(varname, strand, LOCATION, LOCSTART, LOCEND,.keep_all=TRUE) 
}
write_xlsx(assoc_snp_annot, filedate("DArT_Assoc_SNPs", ".xlsx", "output/snippy_multi", dateformat = FALSE), sheet = "DArT_Assoc_SNPs", 
           asTable=FALSE)

assoc_snp_tx <- transcripts(txdb, filter=list("gene_id"=assoc_snp_annot$gene_id),
                            columns=c("tx_name")) %>%   as.data.frame(.) # %>% mutate(varname=paste(seqnames, start, sep="_")) , "gene_id"
fa <- open(FaFile(genome_fa))
cds_seqs <- extractTranscriptSeqs(fa,
                                  cdsBy(txdb, by="tx", use.names=TRUE))
close(fa)
tx_seqs <- cds_seqs[assoc_snp_tx$tx_name]
# Save associated SNPs transcripts to file
seqinr::write.fasta(as.list(translate(tx_seqs)), names = names(tx_seqs),
            file.out = filedate("Assoc_SNP_genes_cds", ".fasta", "output/snippy_multi", dateformat = FALSE))
# Save all transcripts to file
seqinr::write.fasta(as.list(translate(cds_seqs)), names = names(cds_seqs),
                    file.out = "output/snippy_multi/All_SNP_genes_cds.fasta")
# effectorP results
effP_res <- tibble::tribble(
                     ~tx_name,            ~Prediction,   ~Probability,
              "mRNA6901", "Non-effector", 0.968,
              "mRNA6907", "Non-effector", 0.973
              )

# Load BLAST results
# blast_fields <- unlist(strsplit("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore subj_desc subj_sci_name subj_kingdom", split=" "))  # pcov
# annot_table <- read_tsv("data/Assoc_SNPs/Assoc_SNP_genes_cds.nr_blastp.outfmt6",
#                 col_names = blast_fields)
loci_annotation <- readxl::read_excel('data/Assoc_SNPs/Assoc_SNP_table_annotated_05_02_2018.xlsx', 
                                      "Tx_annotation") %>%
  group_by(qseqid) %>% summarise(gene_description=paste(subj_desc, collapse=" | ")) %>%
  left_join(as.data.frame(assoc_snp_tx), ., by=c("tx_name"="qseqid")) %>%
  mutate(gene_id=as.character(gene_id))
# columns(txdb)
assoc_table <- left_join(assoc_snp_annot, loci_annotation[c("gene_id", "gene_description")]) #  %>% 
  # dplyr::select(-name) %>% as_tibble() %>% distinct()       #  by=c("GENEID"="gene_id"))
write_tsv(assoc_table, filedate("Assoc_SNP_table_annotated", ".txt", "data/Assoc_SNPs"))
# xlsx::write.xlsx(as.data.frame(assoc_table), 'data/Assoc_SNPs/Assoc_SNP_table_annotated_05_02_2018.xlsx', 
                 # "SNP_annotation", row.names = FALSE, append = TRUE)
# Save data
save.image(file = filedate("association_analysis_DART", ".RData", "output/snippy_multi/"))

load("data/association_analysis_all_15_02_2018.RData")

# Visualise as a circos
#ctg13 <- genome_seqs$ctg13

# Sum up number of SNPs in bins
# seqlengths <- setNames(circ_data$End, circ_data$Chrom)
# bins <- tileGenome(seqlengths, tilewidth=10000, cut.last.tile.in.chrom=TRUE)
# tileranges <- unlist(tileGenome(seqinfo(allvars), tilewidth=5000))
# hits.df <- as.data.frame(findOverlaps(tileranges, allvars))
# unique(hits.df$queryHits)

##### Circlize plot  #####
# Create data with "chromosome" information
circ_data <- data.frame(Chrom=names(genome_seqs), Start=1,
                        End=sapply(genome_seqs, length)) %>% filter(End>=5e5)
# Convert SNP association data to bed format
assoc_data <- signif_snps[c("chr", "pos", "pvalue", "qvalue")] %>%
  mutate(end=pos, trans_qvalue=-log10(qvalue), trans_pvalue=-log10(pvalue)) %>%   .[c("chr", "pos", "end", "pvalue", "trans_pvalue")] %>% mutate_at(vars(pos, end), .funs = as.numeric) #, "qvalue", "trans_qvalue")]
FDR <- 0.05

# Combine all variants from all samples, then split to a list by LOCATION
# Prepare track data as a list of data frames
bed_list <- varsDF %>% .[2:7]  %>% droplevels(.) %>% split(., f=.$LOCATION)
# Add each Location as a separate track

cols=brewer.pal(length(bed_list), "Set1")

# prepare legends
# discrete
pdf("plots/Circos_legend.pdf", width = 5, height = 5)
lgd_points = ComplexHeatmap::Legend(at = pretty_term(names(bed_list)), type = "points", pch=17,
                    legend_gp = gpar(col = cols, fill="white"), title_position = "topleft",
                    title = "SNP Location")

lgd_lines = ComplexHeatmap::Legend(at = c(sprintf("Significant (p-value<%.3f)",p_val), "Non-signficant"),type = "points", pch=124, #type = "lines",
                   legend_gp = gpar(col = c("darkred", "darkslategrey"), size = 2), title_position = "topleft",
                   title = "SNP Association")
grid.arrange(lgd_points, lgd_lines, ncol=1 )
dev.off()

pdf(filedate("Circos_plot_SNPs_outside", ".pdf", "plots"), width = 8, height=8)
# Initialize circluar plot
circos.par(start.degree=90, gap.degree=c(rep(2, nrow(circ_data)-1), 3)) # increase gap between first and last and add y-axis
circos.genomicInitialize(circ_data)

for (i in 1:length(bed_list)) {
  # determine X and Y to plot text in the first cell

  circos.genomicDensity(bed_list[[i]], window.size = 1e5, col=cols[i], track.height = 0.08)
  i_track = get.cell.meta.data("track.index")
  # circos.yaxis("left", sector.index="ctg00", track.index = i_track, labels = NA)
  # X = get.cell.meta.data("xlim", sector.index = "ctg01", track.index = i_track)
  # Y = get.cell.meta.data("ylim", sector.index = "ctg01", track.index = i_track)
  #
  # circos.text(X[2]*0.5, Y[2], names(bed_list)[i], cex = 0.8, col=cols[i],track.index=i_track,
  #             sector.index="ctg01" ,facing = "bending.inside", niceFacing = TRUE)

}
# Add a track to show SNP association with pathogenicity
circos.genomicTrackPlotRegion(assoc_data, ylim = c(0, 8) , track.height = 0.08,
                              panel.fun = function(region, value, ...) {
                            line_cols <- dplyr::if_else(value$pvalue<=p_val, "darkred", "darkslategrey")
                                  # sapply(value[1], 
                                  #             function(p) ifelse(p<=p_val, "darkred", "darkslategrey"))
                                circos.genomicLines(region, value[2], type = "h", col=line_cols)
                              })
dev.off()
circos.clear()



# i_track = get.cell.meta.data("track.index")
# X = get.cell.meta.data("xlim", sector.index = "ctg01", track.index = i_track)
# Y = get.cell.meta.data("ylim", sector.index = "ctg01", track.index = i_track)
#
# circos.text(X[2]*0.5, Y[2], "q-value", cex = 0.8, col="black",track.index=i_track,
#             sector.index="ctg01" ,facing = "bending.inside", niceFacing = TRUE)


# Produce a circular plot to show SNP association with
circos.genomicTrackPlotRegion(test, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                i=getI(...)
                                circos.genomicPoints(region, value, ..., cex=0.35, pch=15,
                                                     col = my_pal[value$Nucl_int])
                                circos.yaxis(at=i, labels=names(test))
                              })

circos.clear()

circos.yaxis(names(test))
circos.yaxis(at, labels, sector.index, track.index)
