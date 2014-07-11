Detection and filtering of fusions with RNA-seq
====================================================

###Software requirements
+ [defuse](http://sourceforge.net/projects/defuse/)
+ [tophat-fusion](http://ccb.jhu.edu/software/tophat/fusion_index.html)
+ [SOAPfuse](http://soap.genomics.org.cn/soapfuse.html)
+ [Additional tools] (https://docs.google.com/spreadsheet/ccc?key=0ArsHWemp6jw_dGlheGZwT21ONjl0WW9VYVEwWEpyYUE#gid=2)

A review/collection of fusion detection tools/algorithm 
+ [Application of next generation sequencing to human gene fusion detection: computational tools, features and perspectives](http://bib.oxfordjournals.org/content/14/4/506)
+ [The structure of state-of-art gene fusion-finder algorithms](https://www.oapublishinglondon.com/article/617)

Benchmark publications
+ [State-of-the-art fusion-finder algorithms sensitivity and specificity.](http://www.ncbi.nlm.nih.gov/pubmed/23555082)
+ [http://www.ncbi.nlm.nih.gov/pubmed/23815381](http://www.ncbi.nlm.nih.gov/pubmed/23815381)

###Sample collection (inhouse)
[google spreadsheet with rna-seq data of granulocytes](linktospreadsheet)


###Post-filtering steps
+ Filtering by the number of supporting reads
+ Removing known read-through events (use AceView database)
+ Removing pseudogenes
+ Removing ribosomal genes
+ Filtering out fusions which support reads overlap with repetitive elements
+ Removing fusions, which partner genes belong to the same family
+ Evidence from visualization (IGV browser)

Summary table of filters used by fusion detection algorithms
![image](../img/filters.png)
