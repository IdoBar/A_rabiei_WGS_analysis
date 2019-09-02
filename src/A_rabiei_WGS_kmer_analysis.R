pacman::p_load(paletteer, scales)
kmer_hist <- read_tsv("raw_data/khist_summary.txt") %>% 
  set_names(., c("Sample", "kmer", "Depth", "Count", "logScale"))
# limit range (ignore unique kmers with very low depth)
plot_hist <- kmer_hist %>% filter(Depth<200, Depth>2, kmer==31) # 

kmer_peaks <- read_tsv("raw_data/kmer_peaks_summary.txt")
peaks_summary <- kmer_peaks %>% filter(k==31) %>% rename(Sample=sample)

# plot for each kmer
ggplot(plot_hist, mapping = aes(x=Depth, y=Count, colour=Sample)) + scale_y_continuous(labels = comma) +
  geom_line(size=1) + scale_color_paletteer_d(ggsci, default_aaas) + # RColorBrewer, Set1; awtools, mpalette
  facet_wrap(vars(Sample)) + plot_theme() + 
  geom_text(aes(x=25, y=1.75e6, 
                label=sprintf("Predicted ploidy=%sn\nGenome size=%s", ploidy, comma(genome_size))), 
                        # format(genome_size, scientific = TRUE, digits = 3))), 
            data=peaks_summary, hjust=0, colour="black", size=4.5)

ggsave("output/haplo.bt2.fb/plots/A_rabiei_2018_kmer_analysis_k31.pdf", width=12, height=7)
  
  group_by(sample) %>% 
  summarise(mean_genome_size=mean(genome_size),genome_size_SE=se(genome_size))
