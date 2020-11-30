# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# load packages we're going to use
pacman::p_load(tidyverse, vcfR, R.utils, RColorBrewer, ggrepel, glue, here)

# read the VCF file
vcf_file <- "../FB_vars_BT2_07_04_2019/A_rabiei_2018_isolates_haplo.bt2.fb.U20000DP.poly.snps.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
# extract depth
dp <- extract.gt(vcf, element="DP", as.numeric=TRUE)
dp[1:10,1:10]
# tidy our data
tidy_dp <- dp  %>% as.data.frame() %>% 
  rownames_to_column(var="loci") %>%
  pivot_longer(cols = -starts_with("loci"), names_to = "genotype", values_to = "dp") 
# visualise distribution
extrafont::loadfonts()
ggplot(tidy_dp, aes(x=genotype, y=dp)) +
  geom_violin(fill="dodgerblue3", alpha=0.5, trim = TRUE) +
  geom_boxplot(width = 0.2, colour="grey30", outlier.size = 1) +
  plot_theme("bw", 18) +
  labs(y="Depth", x="Genotype") +
  theme_niwot() +
  theme(axis.text.x = element_text(angle = -30, hjust=0, size = rel(0.6)))
ggsave("output/depth_analysis/depth_distribution_per_sample.pdf", width = 15, height = 10)
