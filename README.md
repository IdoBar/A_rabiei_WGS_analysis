-   [Experimental Design](#experimental-design)
-   [Aims](#aims)
-   [Analysis Pipeline](#analysis-pipeline)
    -   [General overview:](#general-overview)
    -   [Methods](#methods)
    -   [Appendix 1. Useful resources](#appendix-1.-useful-resources)
    -   [Appendix 2. General
        information](#appendix-2.-general-information)
    -   [Bibliography](#bibliography)

Experimental Design
===================

In 2017, DNA was extracted from 21 strains of *Ascochyta rabiei* and
sent for Whole-Genome-Sequencing (WGS) on an Illumina HiSeq2500,
producing 100 bp short paired-end reads (Macrogen, Korea).  
In the following year (2018), DNA from 20 additional *A. rabiei*
isolates was extracted and sent for WGS, first to AgriBio, Centre for
AgriBioscience, Agriculture Victoria Research and on a HiSeq2500,
producing 150 bp short paired-end reads. Since the library preparation
and sequencing was substantially delayed, 18 DNA samples, mostly
overlapping with the 20 samples sent for AgriVic, were sent for
sequencing at the Australian Genome Research Facility (AGRF, Melbourne)
on 4 lanes of a NextSeq flowcell, producing 150 bp paried-end reads (run
name CAGRF19461).  
Details of the sequenced isolates is provided in (Table 1).

    datatable(as.data.frame(samples_table), caption=tbls("samples")) %>% # , 
              # options = list(dom = 'tf', pageLength = 40)) %>%
      formatStyle('Pathogenicity',
      backgroundColor = styleInterval(0:3, c('limegreen','gold', 'orange', 'orangered', 'firebrick'))
    )# pander , justify="left"

<!--html_preserve-->

<script type="application/json" data-for="htmlwidget-5b0daba8c863f041a4fc">{"x":{"filter":"none","caption":"<caption>Table  1: Ascochyta rabiei isolates used for DNA sequencing.<\/caption>","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40"],["TR9529","TR9571","TR9573","F17191-1","16CUR018","15CUR002","15CUR005","TR6417","FT13092-2","17CUR007","F17076-2","TR9543","TR9568","16CUR017","16CUR019","F16083-1","F16253-1","15DON007","FT15023","FT15025","FT15028","FT15029","FT15030","FT13092-4","16CUR015","TR8102","15CUR001","FT13092-6","16RUP012","16RUP013","TR8105","15CUR003","14DON003","TR6400","F17067-1","F17175-1","TR9544","TR9538","15DON001","TR6408"],["Chinchilla","Gurley","Gurley","Pt Broughton","Curyo","Curyo","Curyo","Yallaroi","Kingsford","Curyo","Finley","Fox Holes","Gurley","Curyo","Curyo","Moonta","Pt Broughton","Donald","Moonta","Moonta","Weetula","Weetula","Weetula","Kingsford","Curyo","Narromine","Curyo","Kingsford","Rupanyup","Rupanyup","Strathdoon, Narromine","Curyo","Donald","Yallaroi","Coonalpyn","Elmore","Fox Holes","Gravel Pit Hill","Donald","Yallaroi"],["QLD","NSW","NSW","SA","VIC","VIC","VIC","NSW","SA","VIC","NSW","QLD","NSW","VIC","VIC","SA","SA","VIC","SA","SA","SA","SA","SA","SA","VIC","NSW","VIC","SA","VIC","VIC","NSW","VIC","VIC","NSW","SA","VIC","QLD","QLD","VIC","NSW"],[2017,2017,2017,2017,2016,2015,2015,2014,2013,2017,2017,2017,2017,2016,2016,2016,2016,2015,2015,2015,2015,2015,2015,2013,2016,2016,2015,2013,2016,2016,2016,2015,2014,2014,2017,2017,2017,2017,2015,2014],["PBA Seamer","PBA Seamer","PBA Seamer","Genesis090","Genesis090","Genesis090","Genesis090","PBA HatTrick","Genesis090","Genesis090","Genesis090","PBA Seamer","PBA Seamer","Genesis090","Genesis090","Genesis090","Genesis090","Slasher","Genesis090","Genesis090","Genesis090","Genesis090","Genesis090","Genesis090","Genesis090","PBA HatTrick","Genesis090","Genesis090","Genesis090","Genesis090","PBA HatTrick","Genesis090","Slasher","PBA HatTrick","Genesis090","Genesis090","PBA Seamer","PBA Seamer","Genesis090","PBA HatTrick"],["High","High","High","High","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low"],["High","High","High","High","High","High","High","High","High","Moderate","High","High","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Moderate","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low","Low"],["High","High","High","High","High","High","High","High","High","High","Moderate","Moderate","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","Low","Low","Low","Low","Low","Low"],["High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","High","Low","Low","Low","Low","Low","Low"],[4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0],[null,null,null,null,null,"ARH007","ARH180","ARH001","ARH001",null,null,null,null,"ARH001","ARH001",null,null,"ARH001","ARH001","ARH001","ARH077","ARH001","ARH001","ARH001","ARH136","ARH001","ARH001","ARH001","ARH001","ARH001",null,"ARH001","ARH001","ARH001",null,null,null,null,"ARH001","ARH001"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Isolate<\/th>\n      <th>Site<\/th>\n      <th>State<\/th>\n      <th>Collection_Year<\/th>\n      <th>Host Cultivar<\/th>\n      <th>ICC3996<\/th>\n      <th>Genesis090<\/th>\n      <th>HatTrick<\/th>\n      <th>Rating<\/th>\n      <th>Pathogenicity<\/th>\n      <th>Haplotype<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false,"rowCallback":"function(row, data) {\nvar value=data[10]; $(this.api().cell(row, 10).node()).css({'background-color':isNaN(parseFloat(value)) ? '' : value <= 0.000000 ? 'limegreen' : value <= 1.000000 ? 'gold' : value <= 2.000000 ? 'orange' : value <= 3.000000 ? 'orangered' : 'firebrick'});\n}"}},"evals":["options.rowCallback"],"jsHooks":[]}</script>
<!--/html_preserve-->
Aims
====

-   Identify strain-unique variants to develop detection methods
-   Associate aggressiveness with specific variants

Analysis Pipeline
=================

General overview:
-----------------

1.  Data pre-processing:
    1.  Quality check
    2.  Adaptor trimming
    3.  Post-trim quality check
2.  Mapping reads to a reference genome (keep unmapped)
3.  Reads deduplication
4.  Variant calling and filtration
5.  Variant annotation
6.  Variant-Pathogenicity association
7.  Produce variant statistics and assessment

Methods
-------

DNA-Seq data processing and genome assembly were performed following the
methods specified by Hagiwara *et al.* (2014) (see details in Appendix
2), Haas *et al.* (2011), Hittalmani *et al.* (2016) and Verma *et al.*
(2016), with modification to run on the *Griffith University Gowonda HPC
Cluster* (using PBSPro scheduler) and using tools from
<font face='consolas'>BBtools</font> (v38.22; Bushnell (2014)). See
official download page on
[SourceForge](https://sourceforge.net/projects/bbmap/), [user
guide](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
and [SEQanswers
thread](http://seqanswers.com/forums/showthread.php?t=42776).  
Detailed methods are provided in the included
`A_rabiei_genome_sequencing.html` file.

Appendix 1. Useful resources
----------------------------

-   Whole-Genome Comparison of *Aspergillus fumigatus* Strains Serially
    Isolated from Patients with Aspergillosis. (Hagiwara *et al.* 2014):

> **Sequence analysis:** The Illumina data sets were trimmed using
> fastq-mcf in ea-utils (version 1.1.2-484), i.e., sequencing adapters
> and sequences with low quality scores (Phred score \[Q\], &lt;30) were
> removed (24). The data sets were mapped to the genome sequence of the
> *A. fumigatus* genome reference strain Af293 (29,420,142 bp, genome
> version s03-m04-r03) (25, 26) using Bowtie 2 (version 2.0.0-beta7)
> with the very sensitive option in end-to-end mode (27). Duplicated
> reads were removed using Picard (version 1.112)
> (<http://picard.sourceforge.net>). The programs mpileup and bcftools
> from SAMtools (version 0.1.19-44428cd) were used to perform further
> quality controls. In mpileup, the -q20 argument was used to trim reads
> with low-quality mapping, whereas the argument -q30 was used to trim
> low-quality bases at the 3' end (28). The bcftools setting was set to
> -c in order to call variants using Bayesian inference. Consensus and
> single nucleotide polymorphisms (SNPs) were excluded if they did not
> meet a minimum coverage of 5x or if the variant was present in &lt;90%
> of the base calls (29, 30). The genotype field in the variant call
> format (VCF) files indicates homozygote and heterozygote probabilities
> as Phred-scaled likelihoods. SNPs were excluded if they were called as
> heterozygous genotypes using SAMtools. The mapping results were
> visualized in the Integrative Genomics Viewer (version 2.3.3) (31,
> 32). The reference genome data included information on open reading
> frames and annotations, from which the SNPs were designated
> non-synonymous or synonymous.  
> Single nucleotide mutations were confirmed by Sanger sequencing.
> Regions of approximately 400 bp that contained a mutation were
> amplified with appropriately designed primer pairs and then sequenced.
> The primer sequences are listed in Table S1 in the supplemental
> material, which were named as follows. For verification of the SNPs in
> strains from patient I or patient II, PaI or PaII was added to the
> primer name, respectively. For non-synonymous SNPs, synonymous SNPs,
> or SNPs in a non-coding region, (NS, Syno, NonC) was added to the
> primer name, respectively.  
> **Analysis of unmapped reads:** *De novo* assembly of the unmapped
> reads was conducted using the Newbler assembler 2.9 (Roche), with
> default parameters. The contigs were selected based on size/depth
> criteria: those of &lt;500 bp and/or with a depth of &lt;30x coverage
> were removed. To investigate whether unique genome sequences were
> present in strains isolated from the same patient, the unmapped reads
> of each strain were mapped to the contigs generated from all the
> strains in the same patient by the Bowtie 2 software. The coverage of
> the mapped regions was then evaluated. Gene predictions were performed
> using the gene prediction tool AUGUSTUS (version 2.5.5), with a
> training set of *A. fumigatus* (33). The parameters of AUGUSTUS were
> -species = aspergillus\_fumigatus, -strand = both, -genemodel =
> partial, -singlestrand = false, -protein = on, -introns = on, -start =
> on, -stop = on, -cds = on, and -gff3 = on. To compare all the
> predicted genes with *Aspergillus* genes, consisting of 244,811 genes
> available on AspGD (34), a reciprocal BLAST best hit approach was
> performed by BLASTp (35), with an E value of 1.0e<sup>-4</sup>. All
> BLASTp results were filtered based on a BLASTp identity of ≥80% and an
> aligned length coverage of ≥80%.

Appendix 2. General information
-------------------------------

This document was last updated at 2019-02-06 00:40:47 using R Markdown
(built with R version 3.5.1 (2018-07-02)). Markdown is a simple
formatting syntax for authoring HTML, PDF, and MS Word documents. It is
especially powerful at authoring documents and reports which include
code and can execute code and use the results in the output. For more
details on using R Markdown see <http://rmarkdown.rstudio.com> and
[Rmarkdown
cheatsheet](https://www.rstudio.com/wp-content/uploads/2016/03/rmarkdown-cheatsheet-2.0.pdf).

------------------------------------------------------------------------

Bibliography
------------

<!-- ```{r results='asis', eval=TRUE} -->
<!-- PrintBibliography(biblio) -->
<!-- ``` -->
Bushnell B (2014) BBMap: A fast, accurate, splice-aware aligner.

Haas BJ, Zeng Q, Pearson MD, et al. (2011) Approaches to Fungal Genome
Annotation. Mycology 2:118–141. doi:
[10.1080/21501203.2011.606851](https://doi.org/10.1080/21501203.2011.606851)

Hagiwara D, Takahashi H, Watanabe A, et al. (2014) Whole-Genome
Comparison of *Aspergillus* *Fumigatus* Strains Serially Isolated from
Patients with Aspergillosis. J Clin Microbiol 52:4202–4209. doi:
[10.1128/JCM.01105-14](https://doi.org/10.1128/JCM.01105-14)

Hittalmani S, Mahesh HB, Mahadevaiah C, Prasannakumar MK (2016) De novo
genome assembly and annotation of rice sheath rot fungus *Sarocladium*
oryzae reveals genes involved in Helvolic acid and Cerulenin
biosynthesis pathways. BMC Genomics 17:271. doi:
[10.1186/s12864-016-2599-0](https://doi.org/10.1186/s12864-016-2599-0)

Verma S, Gazara RK, Nizam S, et al. (2016) Draft genome sequencing and
secretome analysis of fungal phytopathogen *Ascochyta* *rabiei* provides
insight into the necrotrophic effector repertoire. Sci Rep 6:
