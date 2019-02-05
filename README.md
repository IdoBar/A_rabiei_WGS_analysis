-   [Whole Genome Sequencing of <i>Ascochyta rabiei</i>
    Isolates](#whole-genome-sequencing-of-ascochyta-rabiei-isolates)
    -   [Experimental Design](#experimental-design)
    -   [Aims](#aims)
    -   [Analysis Pipeline](#analysis-pipeline)
        -   [General overview:](#general-overview)
        -   [Methods](#methods)
    -   [Appendices](#appendices)
        -   [Appendix 1. Useful
            resources](#appendix-1.-useful-resources)
        -   [Appendix 2. General
            information](#appendix-2.-general-information)
    -   [Bibliography](#bibliography)

Whole Genome Sequencing of <i>Ascochyta rabiei</i> Isolates
===========================================================

Experimental Design
-------------------

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

    pander(as.data.frame(samples_table), caption=tbls("samples"), justify="left") 

<table>
<caption>Table 1: Ascochyta rabiei isolates used for DNA sequencing. (continued below)</caption>
<colgroup>
<col width="15%" />
<col width="31%" />
<col width="10%" />
<col width="23%" />
<col width="19%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Isolate</th>
<th align="left">Site</th>
<th align="left">State</th>
<th align="left">Collection_Year</th>
<th align="left">Host Cultivar</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">TR9529</td>
<td align="left">Chinchilla</td>
<td align="left">QLD</td>
<td align="left">2017</td>
<td align="left">PBA Seamer</td>
</tr>
<tr class="even">
<td align="left">TR9571</td>
<td align="left">Gurley</td>
<td align="left">NSW</td>
<td align="left">2017</td>
<td align="left">PBA Seamer</td>
</tr>
<tr class="odd">
<td align="left">TR9573</td>
<td align="left">Gurley</td>
<td align="left">NSW</td>
<td align="left">2017</td>
<td align="left">PBA Seamer</td>
</tr>
<tr class="even">
<td align="left">F17191-1</td>
<td align="left">Pt Broughton</td>
<td align="left">SA</td>
<td align="left">2017</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">16CUR018</td>
<td align="left">Curyo</td>
<td align="left">VIC</td>
<td align="left">2016</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">15CUR002</td>
<td align="left">Curyo</td>
<td align="left">VIC</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">15CUR005</td>
<td align="left">Curyo</td>
<td align="left">VIC</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">TR6417</td>
<td align="left">Yallaroi</td>
<td align="left">NSW</td>
<td align="left">2014</td>
<td align="left">PBA HatTrick</td>
</tr>
<tr class="odd">
<td align="left">FT13092-2</td>
<td align="left">Kingsford</td>
<td align="left">SA</td>
<td align="left">2013</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">17CUR007</td>
<td align="left">Curyo</td>
<td align="left">VIC</td>
<td align="left">2017</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">F17076-2</td>
<td align="left">Finley</td>
<td align="left">NSW</td>
<td align="left">2017</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">TR9543</td>
<td align="left">Fox Holes</td>
<td align="left">QLD</td>
<td align="left">2017</td>
<td align="left">PBA Seamer</td>
</tr>
<tr class="odd">
<td align="left">TR9568</td>
<td align="left">Gurley</td>
<td align="left">NSW</td>
<td align="left">2017</td>
<td align="left">PBA Seamer</td>
</tr>
<tr class="even">
<td align="left">16CUR017</td>
<td align="left">Curyo</td>
<td align="left">VIC</td>
<td align="left">2016</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">16CUR019</td>
<td align="left">Curyo</td>
<td align="left">VIC</td>
<td align="left">2016</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">F16083-1</td>
<td align="left">Moonta</td>
<td align="left">SA</td>
<td align="left">2016</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">F16253-1</td>
<td align="left">Pt Broughton</td>
<td align="left">SA</td>
<td align="left">2016</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">15DON007</td>
<td align="left">Donald</td>
<td align="left">VIC</td>
<td align="left">2015</td>
<td align="left">Slasher</td>
</tr>
<tr class="odd">
<td align="left">FT15023</td>
<td align="left">Moonta</td>
<td align="left">SA</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">FT15025</td>
<td align="left">Moonta</td>
<td align="left">SA</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">FT15028</td>
<td align="left">Weetula</td>
<td align="left">SA</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">FT15029</td>
<td align="left">Weetula</td>
<td align="left">SA</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">FT15030</td>
<td align="left">Weetula</td>
<td align="left">SA</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">FT13092-4</td>
<td align="left">Kingsford</td>
<td align="left">SA</td>
<td align="left">2013</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">16CUR015</td>
<td align="left">Curyo</td>
<td align="left">VIC</td>
<td align="left">2016</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">TR8102</td>
<td align="left">Narromine</td>
<td align="left">NSW</td>
<td align="left">2016</td>
<td align="left">PBA HatTrick</td>
</tr>
<tr class="odd">
<td align="left">15CUR001</td>
<td align="left">Curyo</td>
<td align="left">VIC</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">FT13092-6</td>
<td align="left">Kingsford</td>
<td align="left">SA</td>
<td align="left">2013</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">16RUP012</td>
<td align="left">Rupanyup</td>
<td align="left">VIC</td>
<td align="left">2016</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">16RUP013</td>
<td align="left">Rupanyup</td>
<td align="left">VIC</td>
<td align="left">2016</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">TR8105</td>
<td align="left">Strathdoon, Narromine</td>
<td align="left">NSW</td>
<td align="left">2016</td>
<td align="left">PBA HatTrick</td>
</tr>
<tr class="even">
<td align="left">15CUR003</td>
<td align="left">Curyo</td>
<td align="left">VIC</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">14DON003</td>
<td align="left">Donald</td>
<td align="left">VIC</td>
<td align="left">2014</td>
<td align="left">Slasher</td>
</tr>
<tr class="even">
<td align="left">TR6400</td>
<td align="left">Yallaroi</td>
<td align="left">NSW</td>
<td align="left">2014</td>
<td align="left">PBA HatTrick</td>
</tr>
<tr class="odd">
<td align="left">F17067-1</td>
<td align="left">Coonalpyn</td>
<td align="left">SA</td>
<td align="left">2017</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">F17175-1</td>
<td align="left">Elmore</td>
<td align="left">VIC</td>
<td align="left">2017</td>
<td align="left">Genesis090</td>
</tr>
<tr class="odd">
<td align="left">TR9544</td>
<td align="left">Fox Holes</td>
<td align="left">QLD</td>
<td align="left">2017</td>
<td align="left">PBA Seamer</td>
</tr>
<tr class="even">
<td align="left">TR9538</td>
<td align="left">Gravel Pit Hill</td>
<td align="left">QLD</td>
<td align="left">2017</td>
<td align="left">PBA Seamer</td>
</tr>
<tr class="odd">
<td align="left">15DON001</td>
<td align="left">Donald</td>
<td align="left">VIC</td>
<td align="left">2015</td>
<td align="left">Genesis090</td>
</tr>
<tr class="even">
<td align="left">TR6408</td>
<td align="left">Yallaroi</td>
<td align="left">NSW</td>
<td align="left">2014</td>
<td align="left">PBA HatTrick</td>
</tr>
</tbody>
</table>

<table style="width:99%;">
<colgroup>
<col width="15%" />
<col width="18%" />
<col width="15%" />
<col width="12%" />
<col width="22%" />
<col width="15%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">ICC3996</th>
<th align="left">Genesis090</th>
<th align="left">HatTrick</th>
<th align="left">Rating</th>
<th align="left">Pathogenicity</th>
<th align="left">Haplotype</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">4</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">4</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">4</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">4</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">4</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">4</td>
<td align="left">ARH007</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">4</td>
<td align="left">ARH180</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">4</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">4</td>
<td align="left">ARH001</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">ARH001</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">ARH001</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">ARH077</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">ARH001</td>
</tr>
<tr class="even">
<td align="left">Moderate</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">3</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Low</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">2</td>
<td align="left">ARH136</td>
</tr>
<tr class="even">
<td align="left">Low</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">2</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Low</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">2</td>
<td align="left">ARH001</td>
</tr>
<tr class="even">
<td align="left">Low</td>
<td align="left">Moderate</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">2</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">1</td>
<td align="left">ARH001</td>
</tr>
<tr class="even">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">1</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">1</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">1</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">1</td>
<td align="left">ARH001</td>
</tr>
<tr class="even">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">High</td>
<td align="left">High</td>
<td align="left">1</td>
<td align="left">ARH001</td>
</tr>
<tr class="odd">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">0</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">0</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">0</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">0</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">0</td>
<td align="left">ARH001</td>
</tr>
<tr class="even">
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">Low</td>
<td align="left">0</td>
<td align="left">ARH001</td>
</tr>
</tbody>
</table>

Aims
----

-   Identify strain-unique variants to develop detection methods
-   Associate aggressiveness with specific variants

Analysis Pipeline
-----------------

### General overview:

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

### Methods

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

Appendices
----------

### Appendix 1. Useful resources

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

### Appendix 2. General information

This document was last updated at 2019-02-06 00:46:15 using R Markdown
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
