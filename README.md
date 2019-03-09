---
title: "Whole Genome Sequencing of <i>Ascochyta rabiei</i> Isolates"
author: "Ido Bar"
date: "18 July 2017"
always_allow_html: yes
output: 
    bookdown::html_document2:
#      css: "style/style.css"
      toc: true
      toc_depth: 3
#      highlight: pygments
#      number_sections: false
    # html_document:
    #   css: "style/style.css"
    #   toc: true
    #   toc_float: true
    #   toc_depth: 3
    #   highlight: pygments
    #   number_sections: false
    #   code_folding: hide
      keep_md: true
bibliography: style/Fungal_genomes.bib
csl: style/springer-basic-improved-author-date-with-italic-et-al-period.csl
---







# Whole Genome Sequencing of <i>Ascochyta rabiei</i> Isolates
## Experimental Design
In 2017, DNA was extracted from 21 strains of _Ascochyta rabiei_ and sent for Whole-Genome-Sequencing (WGS) on an Illumina HiSeq2500, producing 100 bp short paired-end reads (Macrogen, Korea).  
In the following year (2018), DNA from 20 additional *A. rabiei* isolates was extracted and sent for WGS, first to AgriBio, Centre for AgriBioscience, Agriculture Victoria Research and on a HiSeq3000, producing 150 bp short paired-end reads. Since the library preparation and sequencing was substantially delayed, 18 DNA samples, mostly overlapping with the 20 samples sent for AgriVic, were sent for sequencing at the Australian Genome Research Facility (AGRF, Melbourne) on 4 lanes of a NextSeq500 flowcell, producing 150 bp paired-end reads (run name CAGRF19461).  
Details of the sequenced isolates is provided in (Table  1).


<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:sample_table)Table  1: Ascochyta rabiei isolates used for DNA sequencing.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Isolate </th>
   <th style="text-align:left;"> Site </th>
   <th style="text-align:left;"> State </th>
   <th style="text-align:right;"> Collection_Year </th>
   <th style="text-align:left;"> Host Cultivar </th>
   <th style="text-align:left;"> Host </th>
   <th style="text-align:left;"> Region </th>
   <th style="text-align:left;"> Haplotype </th>
   <th style="text-align:left;"> Pathotype </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> TR9529 </td>
   <td style="text-align:left;"> Chinchilla </td>
   <td style="text-align:left;"> QLD </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> PBA Seamer </td>
   <td style="text-align:left;"> Seamer </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8B2252;">Extreme</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR9571 </td>
   <td style="text-align:left;"> Gurley </td>
   <td style="text-align:left;"> NSW </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> PBA Seamer </td>
   <td style="text-align:left;"> Seamer </td>
   <td style="text-align:left;"> Reg3 </td>
   <td style="text-align:left;"> ARH09 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8B2252;">Extreme</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR9573 </td>
   <td style="text-align:left;"> Gurley </td>
   <td style="text-align:left;"> NSW </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> PBA Seamer </td>
   <td style="text-align:left;"> Seamer </td>
   <td style="text-align:left;"> Reg2 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8B2252;">Extreme</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> F17191-1 </td>
   <td style="text-align:left;"> Pt Broughton </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #8B2252;">Extreme</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> F17076-2 </td>
   <td style="text-align:left;"> Finley </td>
   <td style="text-align:left;"> NSW </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #B22222;">Very High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR9543 </td>
   <td style="text-align:left;"> Fox Holes </td>
   <td style="text-align:left;"> QLD </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> PBA Seamer </td>
   <td style="text-align:left;"> Seamer </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #B22222;">Very High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16CUR018 </td>
   <td style="text-align:left;"> Curyo </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #B22222;">Very High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15CUR002 </td>
   <td style="text-align:left;"> Curyo </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH02 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #B22222;">Very High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15CUR005 </td>
   <td style="text-align:left;"> Curyo </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH20 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #B22222;">Very High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR6417 </td>
   <td style="text-align:left;"> Yallaroi </td>
   <td style="text-align:left;"> NSW </td>
   <td style="text-align:right;"> 2014 </td>
   <td style="text-align:left;"> PBA HatTrick </td>
   <td style="text-align:left;"> HatTrick </td>
   <td style="text-align:left;"> Reg3 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #B22222;">Very High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FT13092-2 </td>
   <td style="text-align:left;"> Kingsford </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2013 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg6 </td>
   <td style="text-align:left;"> ARH04 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #B22222;">Very High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 17CUR007 </td>
   <td style="text-align:left;"> Curyo </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR9568 </td>
   <td style="text-align:left;"> Gurley </td>
   <td style="text-align:left;"> NSW </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> PBA Seamer </td>
   <td style="text-align:left;"> Seamer </td>
   <td style="text-align:left;"> Reg3 </td>
   <td style="text-align:left;"> ARH09 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16CUR017 </td>
   <td style="text-align:left;"> Curyo </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16CUR019 </td>
   <td style="text-align:left;"> Curyo </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> F16083-1 </td>
   <td style="text-align:left;"> Moonta </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg6 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> F16253-1 </td>
   <td style="text-align:left;"> Pt Broughton </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15DON007 </td>
   <td style="text-align:left;"> Donald </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Slasher </td>
   <td style="text-align:left;"> Slasher </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FT15023 </td>
   <td style="text-align:left;"> Moonta </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FT15025 </td>
   <td style="text-align:left;"> Moonta </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FT15028 </td>
   <td style="text-align:left;"> Weetula </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FT15029 </td>
   <td style="text-align:left;"> Weetula </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FT15030 </td>
   <td style="text-align:left;"> Weetula </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg6 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FT13092-4 </td>
   <td style="text-align:left;"> Kingsford </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2013 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg6 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: white;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FF4500;">High</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16CUR015 </td>
   <td style="text-align:left;"> Curyo </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH14 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFA500;">Moderate</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR8102 </td>
   <td style="text-align:left;"> Narromine </td>
   <td style="text-align:left;"> NSW </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> PBA HatTrick </td>
   <td style="text-align:left;"> HatTrick </td>
   <td style="text-align:left;"> Reg4 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFA500;">Moderate</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15CUR001 </td>
   <td style="text-align:left;"> Curyo </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFA500;">Moderate</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FT13092-6 </td>
   <td style="text-align:left;"> Kingsford </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2013 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg6 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFA500;">Moderate</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16RUP012 </td>
   <td style="text-align:left;"> Rupanyup </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFD700;">Medium</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16RUP013 </td>
   <td style="text-align:left;"> Rupanyup </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFD700;">Medium</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR8105 </td>
   <td style="text-align:left;"> Strathdoon, Narromine </td>
   <td style="text-align:left;"> NSW </td>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> PBA HatTrick </td>
   <td style="text-align:left;"> HatTrick </td>
   <td style="text-align:left;"> Reg4 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFD700;">Medium</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15CUR003 </td>
   <td style="text-align:left;"> Curyo </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH20 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFD700;">Medium</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 14DON003 </td>
   <td style="text-align:left;"> Donald </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2014 </td>
   <td style="text-align:left;"> Slasher </td>
   <td style="text-align:left;"> Slasher </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFD700;">Medium</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR6400 </td>
   <td style="text-align:left;"> Yallaroi </td>
   <td style="text-align:left;"> NSW </td>
   <td style="text-align:right;"> 2014 </td>
   <td style="text-align:left;"> PBA HatTrick </td>
   <td style="text-align:left;"> HatTrick </td>
   <td style="text-align:left;"> Reg3 </td>
   <td style="text-align:left;"> ARH04 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #FFD700;">Medium</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> F17067-1 </td>
   <td style="text-align:left;"> Coonalpyn </td>
   <td style="text-align:left;"> SA </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #32CD32;">Low</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> F17175-1 </td>
   <td style="text-align:left;"> Elmore </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #32CD32;">Low</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR9544 </td>
   <td style="text-align:left;"> Fox Holes </td>
   <td style="text-align:left;"> QLD </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> PBA Seamer </td>
   <td style="text-align:left;"> Seamer </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #32CD32;">Low</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR9538 </td>
   <td style="text-align:left;"> Gravel Pit Hill </td>
   <td style="text-align:left;"> QLD </td>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> PBA Seamer </td>
   <td style="text-align:left;"> Seamer </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #32CD32;">Low</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15DON001 </td>
   <td style="text-align:left;"> Donald </td>
   <td style="text-align:left;"> VIC </td>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Genesis090 </td>
   <td style="text-align:left;"> Reg5 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #32CD32;">Low</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TR6408 </td>
   <td style="text-align:left;"> Yallaroi </td>
   <td style="text-align:left;"> NSW </td>
   <td style="text-align:right;"> 2014 </td>
   <td style="text-align:left;"> PBA HatTrick </td>
   <td style="text-align:left;"> HatTrick </td>
   <td style="text-align:left;"> Reg3 </td>
   <td style="text-align:left;"> ARH01 </td>
   <td style="text-align:left;"> <span style="     color: black;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #32CD32;">Low</span> </td>
  </tr>
</tbody>
</table>

## Aims
* Identify strain-unique variants to develop detection methods
* Associate aggressiveness with specific variants

## Analysis Pipeline
### General overview:
1. Data pre-processing:
    a. Quality check
    b. Adaptor trimming
    c. Post-trim quality check
2. Mapping reads to a reference genome (keep unmapped)
3. Reads deduplication
4. Variant calling and filtration
5. Variant annotation
6. Variant-Pathogenicity association
7. Produce variant statistics and assessment 

### Methods
DNA-Seq data processing, mapping and variant calling were performed on the _Griffith University Gowonda HPC Cluster_ (using Torque scheduler), following the methods specified by @hagiwara_whole-genome_2014 (see details in Appendix 2), @haas_approaches_2011, @hittalmani_novo_2016 and @verma_draft_2016, with modification to use <font face='consolas'>FreeBayes</font> v1.2.0 [@garrison_haplotype-based_2012] to assign variant probability scores and call variants.  
An alternative approach was tested, using a complete suite of tools from <font face='consolas'>BBtools</font> v38.22; @bushnell_bbmap:_2014. See official download page on [SourceForge](https://sourceforge.net/projects/bbmap/), [user guide](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) and [SEQanswers thread](http://seqanswers.com/forums/showthread.php?t=42776).   
Detailed methods, including code for running each of the analyses steps are provided in the associated [A_rabiei_WGS_analysis GitHub repository](https://github.com/IdoBar/A_rabiei_WGS_analysis).

## Appendices
### Appendix 1. Useful resources

* Whole-Genome Comparison of _Aspergillus fumigatus_ Strains Serially Isolated from Patients with Aspergillosis. [@hagiwara_whole-genome_2014]:

> **Sequence analysis:** The Illumina data sets were trimmed using fastq-mcf in ea-utils (version 1.1.2-484), i.e., sequencing adapters and sequences with low quality scores (Phred score [Q], <30) were removed (24). The data sets were mapped to the genome sequence of the _A. fumigatus_ genome reference strain Af293 (29,420,142 bp, genome version s03-m04-r03) (25, 26) using Bowtie 2 (version 2.0.0-beta7) with the very sensitive option in end-to-end mode (27). Duplicated reads were removed using Picard (version 1.112) (<http://picard.sourceforge.net>). The programs mpileup and bcftools from SAMtools (version 0.1.19-44428cd) were used to perform further quality controls. In mpileup, the -q20 argument was used to trim reads with low-quality mapping, whereas the argument -q30 was used to trim low-quality bases at the 3' end (28). The bcftools setting was set to -c in order to call variants using Bayesian inference. Consensus and single nucleotide polymorphisms (SNPs) were excluded if they did not meet a minimum coverage of 5x or if the variant was present in <90% of the base calls (29, 30). The genotype field in the variant call format (VCF) files indicates homozygote and heterozygote probabilities as Phred-scaled likelihoods. SNPs were excluded if they were called as heterozygous genotypes using SAMtools. The mapping results were visualized in the Integrative Genomics Viewer (version 2.3.3) (31, 32). The reference genome data included information on open reading frames and annotations, from which the SNPs were designated non-synonymous or synonymous.  
Single nucleotide mutations were confirmed by Sanger sequencing. Regions of approximately 400 bp that contained a mutation were amplified with appropriately designed primer pairs and then sequenced. The primer sequences are listed in Table S1 in the supplemental material, which were named as follows. For verification of the SNPs in strains from patient I or patient II, PaI or PaII was added to the primer name, respectively. For non-synonymous SNPs, synonymous SNPs, or SNPs in a non-coding region, (NS, Syno, NonC) was added to the primer name, respectively.  
**Analysis of unmapped reads:** _De novo_ assembly of the unmapped reads was conducted using the Newbler assembler 2.9 (Roche), with default parameters. The contigs were selected based on size/depth criteria: those of <500 bp and/or with a depth of <30x coverage were removed. To investigate whether unique genome sequences were present in strains isolated from the same patient, the unmapped reads of each strain were mapped to the contigs generated from all the strains in the same patient by the Bowtie 2 software. The coverage of the mapped regions was then evaluated. Gene predictions were performed using the gene prediction tool AUGUSTUS (version 2.5.5), with a training set of  _A. fumigatus_ (33). The parameters of AUGUSTUS were -species = aspergillus_fumigatus, -strand = both, -genemodel = partial, -singlestrand = false, -protein = on, -introns = on, -start = on, -stop = on, -cds = on, and -gff3 = on. To compare all the predicted genes with _Aspergillus_ genes, consisting of 244,811 genes available on AspGD (34), a reciprocal BLAST best hit approach was performed by BLASTp (35), with an E value of 1.0e<sup>-4</sup>. All BLASTp results were filtered based on a BLASTp identity of $\ge80$% and an aligned length coverage of $\ge80$%.


### Appendix 2. General information
This document was last updated at 2019-03-10 02:28:05 using R Markdown (built with R version 3.5.1 (2018-07-02)). Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. It is especially powerful at authoring documents and reports which include code and can execute code and use the results in the output. For more details on using R Markdown see <http://rmarkdown.rstudio.com> and [Rmarkdown cheatsheet](https://www.rstudio.com/wp-content/uploads/2016/03/rmarkdown-cheatsheet-2.0.pdf).

***
## Bibliography

<!-- ```{r results='asis', eval=TRUE} -->
<!-- PrintBibliography(biblio) -->
<!-- ``` -->

