rnaseqlib 
=========

A collection of Python modules & Rscripts for analysis of RNA-seq data for the following publication:
### Mutational Landscape of the Transcriptome Offers a Rich Neoantigen Resource for Immunotherapy of Myeloproliferative Neoplasms
Fiorella Schischlik et al. 2018

### Software requirements
#### `Fusion detection`
+ [defuse](http://sourceforge.net/projects/defuse/)
+ [tophat-fusion](http://ccb.jhu.edu/software/tophat/fusion_index.html)
+ [SOAPfuse](http://soap.genomics.org.cn/soapfuse.html)

#### `Variant calling`
+ [STAR](https://github.com/alexdobin/STAR)
+ [Samtools](http://samtools.sourceforge.net/)
+ [GATK](https://www.broadinstitute.org/gatk/)
+ [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
+ [Bedtools](https://github.com/arq5x/bedtools2)

#### `Aberrant splicing`
+ [STAR](https://github.com/alexdobin/STAR)
+ [deboever-sf3b1-2015](https://github.com/cdeboever3/deboever-sf3b1-2015)

#### `External Python modules`
+ [pandas](https://github.com/pydata/pandas) [required]
+ [PyVCF](https://pypi.python.org/pypi/PyVCF) [optional; for variant calling]
+ [gffutils](https://pypi.python.org/pypi/gffutils) [optional; for annotating fusions]
+ [pyBigWig](https://pypi.python.org/pypi/pyBigWig) [optional; for annotating fusions]

#### `Environment variables used`
+ ```$TMPDIR``` (path to temporary directory)  
+ ```$NGS_GATK``` (path to gatk executables) 
+ ```$NGS_PICARD``` (path to picard executables)

Specific workflows exits for:
+ [Annotation and filtering of fusions detected with RNA-seq.](doc/fusion_workflow.md)
+ [Calling and annotating variants with RNA-seq.](doc/variant_calling.md)


## License

rnaseqlib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

rnaseqlib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

