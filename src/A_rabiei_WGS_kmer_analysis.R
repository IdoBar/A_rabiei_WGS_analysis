# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# load packages we're going to use
pacman::p_load(tidyverse, paletteer, scales)
palettes_d_names %>% filter(palette=="Bold" )
palettes_d_names %>% filter(length>8, type=="qualitative")
kmer_hist <- read_tsv("raw_data/khist_summary.txt") %>% 
  set_names(., c("Sample", "kmer", "Depth", "Count", "logScale"))
# kmer_hist %>% group_by(Sample) %>% slice(max(main_peak))

kmer_peaks <- read_tsv("raw_data/kmer_peaks_summary.txt")
# summarise peaks to calculate average predicted genome size and 
peaks_summary <- kmer_peaks %>% group_by(sample) %>% rename(Sample=sample) %>% 
  summarise(mean_genome_size=mean(genome_size),genome_size_SE=se(genome_size), max_cov=max(main_peak), 
            ploidy=max(ploidy)) %>% 
  arrange(max_cov) %>% mutate(Sample=forcats::fct_inorder(Sample), 
                              genome_label=sprintf("%s\u2009bp (±%s)", comma(mean_genome_size/1000), 
                                           round(genome_size_SE, digits = 0)))
# limit range (ignore unique kmers with very low depth)
plot_hist <- kmer_hist %>% filter(Depth<200, Depth>2) %>% 
  mutate(kmer=factor(kmer), Sample=factor(Sample, levels = levels(peaks_summary$Sample)))#, kmer==31) #  

# filter(k==31) %>% rename(Sample=sample)

# plot for each kmer (check out the following to support UniCODE characters: https://stackoverflow.com/questions/12768176/unicode-characters-in-ggplot2-pdf-output) , thin space: \u2009
ggplot(plot_hist, mapping = aes(x=Depth, y=Count, colour=kmer)) + 
  scale_y_continuous(labels = comma) +
  geom_line(size=1) + scale_color_paletteer_d("pals::tol") + # ggsci, default_aaas; RColorBrewer, Set1; awtools, mpalette
  facet_wrap(vars(Sample)) + plot_theme() + 
  geom_text(aes(x=15, y=2.0e6, 
                label=sprintf("Predicted ploidy=%sn\nGenome size=%s(±%s) Kbp", ploidy, 
                              comma(mean_genome_size/1000), round(genome_size_SE/1000, digits = 1))),  
                        # format(genome_size, scientific = TRUE, digits = 3))) / comma(genome_size)
            data=peaks_summary, hjust=0, colour="black", size=5)

ggsave("output/haplo.bt2.fb/plots/A_rabiei_2018_kmer_analysis_size_CI.pdf", width=13, height=7.5)
  
  
