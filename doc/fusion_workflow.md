Detection and filtering of fusions with RNA-seq
====================================================

###Software requirements
+ [defuse](http://sourceforge.net/projects/defuse/)
+ [tophat-fusion](http://ccb.jhu.edu/software/tophat/fusion_index.html)
+ [SOAPfuse](http://soap.genomics.org.cn/soapfuse.html)
+ [Summary of fusion detection tools (google spreadsheet](https://docs.google.com/spreadsheet/ccc?key=0ArsHWemp6jw_dGlheGZwT21ONjl0WW9VYVEwWEpyYUE#gid=2)

A review/collection of fusion detection tools/algorithm 
+ [Application of next generation sequencing to human gene fusion detection: computational tools, features and perspectives](http://bib.oxfordjournals.org/content/14/4/506)
+ [The structure of state-of-art gene fusion-finder algorithms](https://www.oapublishinglondon.com/article/617)

Benchmark publications
+ [State-of-the-art fusion-finder algorithms sensitivity and specificity.](http://www.ncbi.nlm.nih.gov/pubmed/23555082)
+ [State of art fusion-finder algorithms are suitable to detect transcription-induced chimeras in normal tissues?](http://www.ncbi.nlm.nih.gov/pubmed/23815381)

###Sample collection (inhouse)
[google spreadsheet with rna-seq data of granulocytes](linktospreadsheet)


###Post-filtering steps

After read mapping and nominating potential fusion candidates a set of filters is applied based on biological and technical indications.

Scoring filters
+ Filter by the number of encompassing and spanning reads (aka supporting reads)

Technical artefact filters
+ Filtering out fusions which support reads overlap with repetitive elements
+ Removing fusions, which partner genes belong to the same family
+ Evidence from visualization (IGV browser)

Biological filters
+ Removing known read-through events (use AceView database)
+ Removing pseudogenes
+ Removing ribosomal genes


Additional filters

Summary table of filters used by fusion detection algorithms.
+ [The structure of state-of-art gene fusion-finder algorithms. M.Becutti, Aug, 2013](https://www.oapublishinglondon.com/article/617)
![image](../img/filters.png)
