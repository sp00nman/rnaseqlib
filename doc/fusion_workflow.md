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

| DATABASE                  | DESCRIPTION                 | DATE   | COMMENT |
| :------------------------ |:----------------------------|:-------|:--------|
| healthy defuse | n = ? healthy granulocytes processed with defuse | ? | ? |
| healthy tophat-fusion | n = ? healthy granulocytes processed with tophat-fusion | ? | ? | 
| healthy soapfuse | n = ? healthy granulocytes processed with soapfuse | ? | ? | 
| [AceView](http://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/index.html?human) database)|  readthrough  |   ?     | ? |
| [ConjoinG database](http://metasystems.riken.jp/conjoing/) | readthrough |  ? | ? |
| [CACG conjoined genes database](http://cgc.kribb.re.kr/map/) | readthrough | ? | ? | 
| NCBI | readthrough | ? |  ? | ? |
| pseudogene - [GENCODE] (http://www.gencodegenes.org/releases/19.html) | ? | ? | ? | 
| non protein coding - [GENCODE] (http://www.gencodegenes.org/releases/19.html) | ? | ? | ? | 

**Annotation**

| DATABASE                  | DESCRIPTION                 | DATE   | COMMENT |
| :------------------------ |:----------------------------|:-------|:--------|
| [TICdb database](http://www.unav.es/genetica/TICdb/) |    |   ?     | ? |
| [Mitelman Database of Chromosome Aberrations in Cancer](http://cgap.nci.nih.gov/Chromosomes/Mitelman.) |    |   ?     | ? |
| [COSMIC database](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/) |    |   ?     | ? |
| [ChimerDB 2.0 database literature-based annotation](http://ercsb.ewha.ac.kr/FusionGene/) | ? | ? | ? | 


**Additional filters**
+ [Viruses/bacteria/phages genomes database (from the NCBI database)](ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/ (required)

Summary table of filters used by fusion detection algorithms.
+ [The structure of state-of-art gene fusion-finder algorithms. M.Becutti, Aug, 2013](https://www.oapublishinglondon.com/article/617)
![image](../img/filters.png)

###What can be improved?

####PROPOSAL 1
#####Motivation/Literature research
+ "Another important group of gene fusions was associated with breakpoints of low-level copy number changes, involving both gains and deletions. These are interesting in the sense that they represent the types of fusion events leading to gene activation with no association with gene amplifications...Fifth, in the vast majority of the fusions (82%), at least one partner gene was located at a copy number breakpoint as revealed by aCGH, indicating that fusion gene formation is closely associated with unbalanced genomic rearrangements, particularly high-level amplifications" [Identification of fusion genes in breast cancer by paired-end RNA-sequencing.](http://genomebiology.com/content/12/1/R6)

+ [Identification of somatically acquired rearrangements in cancer using genome-wide massively parallel paired-end sequencing.](http://www.nature.com/ng/journal/v40/n6/fig_tab/ng.128_F3.html)


#####Idea
+ Can we use public data/database on germline structural variant calls to filter for gene fusions
  + Most biological database filter for gene IDs (only some databases offer breakpoint positions as well, but most don't)
  + SNP arrays have far too low resolution, CNV information from WGS data could potentially be beneficial for this approach

####PROPOSAL 2
+ Detection and filtering of fusions is followed by manually looking at the sequence and blasting the fusion junction to the reference genome. Very often these sequences fall within repetitive sequences, or the 5' and 3' region of the fusion is very similar

#####Idea
  + Implement a repetitive sequence filter
    + UCSC Mappability Track, Repeat Masker Tracks 
    + [Duke Exclude Regions](http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/) 
    + [DAC Blacklisted Regions](http://hgwdev.cse.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability)
    + [1000 Genomes Masks](http://www.1000genomes.org/announcements/genome-accessibility-information-now-available-1000-genomes-browser-2012-09-06)
  + Implement a 5' 3' homology filter
  + [sequence similarity search](http://www.ebi.ac.uk/Tools/sss/)


