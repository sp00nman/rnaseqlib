Annotation and filtering of fusions with RNA-seq
====================================================

### Annotating fusions
The following python script was designed to annotate fusions and in later steps 
to enrich for novel and somatic fusions. Databases used for annotation are listed
below.


### Usage 
```
usage: annotate_fusions.py [-h] [--project_name PROJECT_NAME]
                           [--input_file INPUT_FILE]
                           [--annotation_file ANNOTATION_FILE]
                           [--output_dir OUTPUT_DIR] [--gtf_dbfile GTF_DBFILE]
                           [--id_conversion ID_CONVERSION]
                           [--min_num_split_reads MIN_NUM_SPLIT_READS]
                           [--min_num_span_reads MIN_NUM_SPAN_READS]
                           [--score SCORE] [--num_cpus NUM_CPUS]

Annotate fusions from fusion tools for RNA-seq.

optional arguments:
  -h, --help            show this help message and exit
  --project_name PROJECT_NAME
                        Name of the project directory.
  --input_file INPUT_FILE
                        Input file with the format (tab separated):
                        tool[defuse,tophatfusion,soapfuse] path/to/file
  --annotation_file ANNOTATION_FILE
                        Annotation file with the format (tab separated):
                        database[see github] /path/to/annotation_file
                        [pairs(A_B)/individual(A,B) [annotate/filter]
  --output_dir OUTPUT_DIR
                        Path to output directory.
  --gtf_dbfile GTF_DBFILE
                        gtf file (transformed to database with gffutils)
  --id_conversion ID_CONVERSION
                        Conversion table. Ensembl ids to genesymbols
  --min_num_split_reads MIN_NUM_SPLIT_READS
                        Minimum amount of split reads.
  --min_num_span_reads MIN_NUM_SPAN_READS
                        Minimum amount of spanning reads.
  --score SCORE         Score (p-value or any other score.)
```

INPUT_FILE (tab separated)

| UNIQ_SAMPLE_ID | TOOL                               | PATH_TO_FUSION_OUTPUTFILE | 
|:----------------|:------------------------------------|:---------------------------|
| eg. MPN0001    | OR [tophatfusion, defuse, soapfuse] | eg. /path/to/file             |

ANNOTATION_FILE (tab separated)

| label | gene_position                               | mode | source | filter_annotate | file_location |
|:----------------|:------------------------------------|:---------------------------|:------|:-----|:-----|
| eg. readthrough    | OR [gene,position] | OR [pair,single,biotype,blacklisted_regions] | eg.conjoing | OR [filter,annotate] | eg. pat/to/annotation_database

GTF_DBFILE

GTF file generated as explained [here](https://pythonhosted.org/gffutils/#create-the-database)

ID_CONVERSION

eg. ENSG00000186716 BCR (tab separated)

--> structure of the database file for source

**pair** eg. ENSG00000186716 ENSG00000143322 (tab separated) 

**single** eg. ENSG00000186716

**biotype** database file created as explained [here](https://pythonhosted.org/gffutils/#create-the-database)

**blacklisted_regions** *.bigWig files


### Annotation databases

| DATABASE                  | DESCRIPTION                 | DATE (last access)   | SOURCE |
| :------------------------ |:----------------------------|:-------|:--------|
| healthy defuse | n=25 healthy individuals; RNA-seq on granulocytes processed with defuse | 2015 | pair |
| healthy tophatfusion | n=25 healthy individuals; RNA-seq on granulocyted processed with tophatfusion | 2015 | pair | 
| healthy soapfuse | n=25 healthy individuals; RNA-seq on granulocytes processed with soapfuse | 2015 | pair | 
| [ConjoinG database](http://metasystems.riken.jp/conjoing/download/) | ConjoinG is a database of 800 Conjoined Genes identified in the Human Genome. |  2009-10-13 | pair |
| [TICdb database](http://www.unav.es/genetica/allseqs_TICdb.txt) | TICdb is a database of Translocation breakpoints In Cancer. It contains 1,374 fusion sequences found in human tumors, involving 431 different genes. |   2013-08     | pair |
| [NCBI](ftp://ftp.ncbi.nih.gov/gene/DATA/gene_group.gz) | Report of genes and their relationships to other genes | 2014_07 | pair |
| 5' or 3' prime biotype gene - [GENCODE] (http://www.gencodegenes.org/releases/19.html) | extract biotype for 5' or 3' gene with GENCODE release 19 | 2015 | biotype | 
| [encode](http://hgwdev.cse.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability) | mapability (50-kmer) | 2016-02 | blacklisted_regions |
| [encode_dac](http://hgwdev.cse.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability) | The DAC Blacklisted Regions aim to identify a comprehensive set of regions in the human genome that have anomalous, unstructured, high signal/read counts in next gen sequencing experiments independent of cell line and type of experiment | 2016-02 | blacklisted regions |
| [encode_duke](http://hgwdev.cse.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability) | The Duke Excluded Regions track displays genomic regions for which mapped sequence tags were filtered out before signal generation and peak calling for Open Chromatin: DNaseI HS and FAIRE tracks| 2016-02 | blacklisted regions |
| [ChimerDB 2.0 database literature-based annotation](http://ercsb.ewha.ac.kr/FusionGene/) | ChimerDB is designed to be a knowledgebase of fusion transcripts collected from various public resources such as the Sanger CGP, OMIM, PubMed, and Mitelman’s database | 2016-02 | pair |
| [HLA] (http://hla.alleles.org/genes/) | List of genes within the HLA Region | 2016-2 | single |
| [HB] (http://ensemble.org) | List of genes part of the hemoglobin complex. | 2016-02 | single | 
| paralogs, fully, partially and same strand overlapping genes | List of genes processed by [fusioncatcher](https://github.com/ndaniel/fusioncatcher) ; ensembl version 72 | 2016-2 | pair |

TODOs:

+ [Mitelman Database of Chromosome Aberrations in Cancer](http://cgap.nci.nih.gov/Chromosomes/Mitelman.)
+ [COSMIC database](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/)
+ [AceView](http://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/index.html?human) database)
+ [CACG conjoined genes database](http://cgc.kribb.re.kr/map/)
+ [Viruses/bacteria/phages genomes database (from the NCBI database)](ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/ (required))
+ [sequence similarity search](http://www.ebi.ac.uk/Tools/sss/)

### Filtering fusions

Fusions are tagged as follows:

| TAG                | DESCRIPTION (fusion is tagged if...) | 
| :------------------|:----------------------------|
| score | defuse score < 0.8 or tophat score < 0 |
| readthrough | fusions appears in any of the above mentioned databases for readthroughs  |
| false_pos | fusions is present in any of these [paralogs, fully, partially and same strand overlapping genes] databases provided by fusioncatcher  |
| healthy | fusion is present in in-house generated databases for fusions detected in healthy individuals|
| no_protein | both fusion partner genes have a biotype marked as non-protein|
| biotype_pseudogene | one of the fusion partner is annotated as pseudogene and has a mapability score of less than 0.5|
| mapability | one of the fusion genes has a mapability score of less then 0.1 |
| blacklisted_region | fusions falls within blacklisted regions |
| distance_lt1000 | fusion partner genes are less 1000 base pairs apart |
| HLA_gene | one of the fusion partners is part of HLA gene family |
| HB_gene | one of the fusion partners is part of HB gene family |
| PASS | fusions passes all of the above filters |
