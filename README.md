rnaseqlib 
=========

A simple collection of modules for analysis of RNA-seq data with focus on:

+ [Annotation and filtering of fusions detected with RNA-seq.](doc/fusion_workflow.md)
+ [Calling and annotating variants with RNA-seq.](doc/variant_calling.md)
+ [Estimating an X-inactivation score (XCI-score) from variants called on X-chromosome with RNA-seq] (doc/xci_skew.mkd)

### Software requirements
--> specific for fusion detection
+ [defuse](http://sourceforge.net/projects/defuse/)
+ [tophat-fusion](http://ccb.jhu.edu/software/tophat/fusion_index.html)
+ [SOAPfuse](http://soap.genomics.org.cn/soapfuse.html)

--> specific for variant calling
+ [STAR] (https://github.com/alexdobin/STAR)
+ [Samtools] (http://samtools.sourceforge.net/)
+ [GATK] (https://www.broadinstitute.org/gatk/)
+ [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
+ [Bedtools](https://github.com/arq5x/bedtools2)

##### Python modules
+ [pandas](https://github.com/pydata/pandas) [required]
+ [PyVCF](https://pypi.python.org/pypi/PyVCF) [optional; for variant calling]
+ [gffutils](https://pypi.python.org/pypi/gffutils) [optional; for annotating fusions]
+ [pyBigWig] (https://pypi.python.org/pypi/pyBigWig) [optional; for annotating fusions]

##### Environment variables
+ ```$TMPDIR``` (path to temporary directory)  
+ ```$NGS_GATK``` (path to gatk executables) 
+ ```$NGS_PICARD``` (path to picard executables)


## License

rnaseqlib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

rnaseqlib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

