Detection and filtering of fusions with RNA-seq
====================================================

###Literature research
A review/collection of fusion detection tools/algorithm 
+ [Application of next generation sequencing to human gene fusion detection: computational tools, features and perspectives](http://bib.oxfordjournals.org/content/14/4/506)
+ [The structure of state-of-art gene fusion-finder algorithms](https://www.oapublishinglondon.com/article/617)

Benchmark publications
+ [State-of-the-art fusion-finder algorithms sensitivity and specificity.](http://www.ncbi.nlm.nih.gov/pubmed/23555082)
+ [State of art fusion-finder algorithms are suitable to detect transcription-induced chimeras in normal tissues?](http://www.ncbi.nlm.nih.gov/pubmed/23815381)

###Software requirements
+ [defuse](http://sourceforge.net/projects/defuse/)
+ [tophat-fusion](http://ccb.jhu.edu/software/tophat/fusion_index.html)
+ [SOAPfuse](http://soap.genomics.org.cn/soapfuse.html)
+ [Summary of fusion detection tools (google spreadsheet](https://docs.google.com/spreadsheet/ccc?key=0ArsHWemp6jw_dGlheGZwT21ONjl0WW9VYVEwWEpyYUE#gid=2)

###Sample collection (inhouse)
[google spreadsheet with rna-seq data of granulocytes](linktospreadsheet)

###Post-filtering steps

After read mapping and nominating potential fusion candidates a set of filters is applied based on biological and technical indications. [Fusioncatcher](https://code.google.com/p/fusioncatcher/wiki/Manual#3.3_-_Genomic_Databases) has a fine collection of databases to filter for.

**Scoring filters**
+ Filter by the number of encompassing and spanning reads (aka supporting reads)

**Technical artefact filters**
+ Filtering out fusions which support reads overlap with repetitive elements
+ Removing fusions, which partner genes belong to the same family
+ Evidence from visualization (IGV browser)

**Biological filters**
+ Removing fusions from healthy individuals
  + [in house](https://docs.google.com/spreadsheet/ccc?key=0ArsHWemp6jw_dGlheGZwT21ONjl0WW9VYVEwWEpyYUE#gid=9) 
+ Removing known/or estimated read-through events
  + [AceView](http://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/index.html?human) database)
  + [ConjoinG database](http://metasystems.riken.jp/conjoing/)
  + [CACG conjoined genes database](http://cgc.kribb.re.kr/map/)
+ Removing pseudogenes
+ Removing ribosomal genes
+ Removing IMGT/HLA genes
  + [IMGT/HLA database](http://www.ebi.ac.uk/ipd/imgt/hla/)
+ Annotating known fusions
  + [COSMIC database](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/)
  + [TICdb database](http://www.unav.es/genetica/TICdb/)
  + [ChimerDB 2.0 database literature-based annotation](http://ercsb.ewha.ac.kr/FusionGene/)
  + [Cancer Genome Project (CGP) translocations database](http://www.sanger.ac.uk/genetics/CGP/Census/)
  + [Mitelman Database of Chromosome Aberrations in Cancer](http://cgap.nci.nih.gov/Chromosomes/Mitelman.)

**Additional filters**
+ [Viruses/bacteria/phages genomes database (from the NCBI database)](ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/ (required)

Summary table of filters used by fusion detection algorithms.
+ [The structure of state-of-art gene fusion-finder algorithms. M.Becutti, Aug, 2013](https://www.oapublishinglondon.com/article/617)
![image](../img/filters.png)

###What can be improved?

#####Literature research
"Another important group of gene fusions was associated with breakpoints of low-level copy number changes, involving both gains and deletions. These are interesting in the sense that they represent the types of fusion events leading to gene activation with no association with gene amplifications...Fifth, in the vast majority of the fusions (82%), at least one partner gene was located at a copy number breakpoint as revealed by aCGH, indicating that fusion gene formation is closely associated with unbalanced genomic rearrangements, particularly high-level amplifications" [Identification of fusion genes in breast cancer by paired-end RNA-sequencing.](http://genomebiology.com/content/12/1/R6)

[Identification of somatically acquired rearrangements in cancer using genome-wide massively parallel paired-end sequencing.](http://www.nature.com/ng/journal/v40/n6/fig_tab/ng.128_F3.html)

[My own data](https://docs.google.com/spreadsheets/d/16tYeYdyrk-Qp4HifYMrcTaeaNtB6xcsLNkbW52Lw-6c/edit#gid=0)

+ Most biological database filter for gene IDs (only some databases offer breakpoint positions as well, but most don't)
  + Can we use public data/database on germline structural variant calls to filter for gene fusions
  + SNP arrays have far too low resolution, CNV information from WGS data could potentially be beneficial for this approach

+ 


