### Workflow for calling variants optimized for RNA sequencing data.

Workflow is well described [here] (http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq)

## Software requirements
 
+ [STAR] (https://github.com/alexdobin/STAR)
+ [Samtools] (http://samtools.sourceforge.net/)
+ [GATK] (https://www.broadinstitute.org/gatk/)
+ [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
+ [Bedtools](https://github.com/arq5x/bedtools2)

##### Python modules
+ [PyVCF]()
+ [pandas](https://github.com/pydata/pandas)

##### Environment variables
+ $TMPDIR (path to temporary directory)  
+ $NGS_GATK (path to gatk executables) 
+ $NGS_PICARD (path to picard executables)

## Usage

```bash
user@node:/path/to/script/rnaseq_varcall.py -h

usage: rnaseq_varcall.py [-h] [--debug DEBUG] [--stage STAGE]
                         [--project_name PROJECT_NAME] [--read1 READ1]
                         [--read2 READ2] [--sample_dir SAMPLE_DIR]
                         [--output_dir OUTPUT_DIR] [--star_genome STAR_GENOME]
                         [--ref_genome REF_GENOME] [--star2pass]
                         [--sample_file SAMPLE_FILE] [--region REGION]
                         [--num_cpus NUM_CPUS] [--annovar ANNOVAR]

Genetic screen workflow 0.0.1

optional arguments:
  -h, --help            show this help message and exit
  --debug DEBUG         Debug level
  --stage STAGE         Limit job submission to a particular analysis stage.[a
                        ll,alignment,extract,replace_rg,duplicates,splitntrim,bqsr,
                        bamfo,samtools,gatk,gatk_filter,hrun,annovar,selection]
  --project_name PROJECT_NAME
                        name of the project
  --read1 READ1         For paired alignment, forward read.
  --read2 READ2         For paired alignment, reverse read.
  --sample_dir SAMPLE_DIR
                        Path to sample directory
  --output_dir OUTPUT_DIR
                        Path to output directory.
  --star_genome STAR_GENOME
                        Genome directory of star aligner.
  --ref_genome REF_GENOME
                        reference genome
  --star2pass           If set, STAR 2-pass mapping will be performed.
  --sample_file SAMPLE_FILE
                        For region variant calling.
  --region REGION       region eg. 20:30946147-31027122
  --num_cpus NUM_CPUS   Number of cpus.
  --annovar ANNOVAR     Annotate variant with annovar.

```

## --stage

### [all]
Process all analysis steps.

### [alignment]
Alignment with STAR and the following additional options set:
```bash
--outFilterIntronMotifs RemoveNoncanonical
--outSAMtype BAM SortedByCoordinate

```

### [star2pass]
If ```star2pass``` is set then a new index is created using splice junction information contained in the file SJ.out.tab from the first pass.

### [replace_rg]
Adds readgroup information and sorts by coordinate as described [here] (http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq)

### [duplicates]
Remove duplicate reads as described [here] (http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq)

### [index]
Index BAM file with samtools.

### [splitntrim]
GATK tool ```-T SplitNCigarReads``` splits reads into exon segments and hard-clip sequences overhanging into the
intronic regions. Additionally mapping qualities are reassigned.

### [indel]
Indel realignment (same as DNA). Recommended [known sites](https://www.broadinstitute.org/gatk/guide/article?id=1247).

### [bqsr]
Base quality recalibration (same as DNA). Recommended [known sites](https://www.broadinstitute.org/gatk/guide/article?id=1247).

### [gatk]
GATK tool ```-T HaplotypeCaller``` caller with options.
```bash
-dontUseSoftClippedBases
-recoverDanglingHeads
-stand_call_conf 20.0

```
Minimum phred-scaled confidence is lowered to 20 according to GATK best practices recommendation.

## Variant Filtering Steps

### [gatk_filter]
GATK ```-T VariantFiltration``` is used to filter for:
- [1] at least 3 SNPs that are within a window of 35 bases ```-window 35 -cluster 3```
- [2] fisher strand value FS>30
- [3] Quality by Depth QD<2

### [hrun]
- [4] filter for variants within homopolymer runs >=5 ```[##FILTER=<ID=HRun,Description="HRun >= 5">]```

### [indel_prox] [not implemented yet]
- [5] filter for SNVs located within 5 bases from an indel (Reumers et al., Nature Biotech., 2012)

### [annovar]
Databases used for annotation [ANNOVAR] (http://annovar.openbioinformatics.org/en/latest/user-guide/download/) and
[inhouse] filtering steps:

| NAME                   | DESCRIPTION                 | DATE   |
| :--------------------- |:----------------------------|:-------|
| cytoBand               | cytoBand annotation             |   ?     |
| phastConsElements46way | ?         |     ?   |
| avsnp142               | dbSNP142 with allelic splitting and left-normalization	| 20141228 |
| snp138NonFlagged       | dbSNP with ANNOVAR index files, after removing those flagged SNPs (SNPs < 1% minor allele frequency (MAF) (or unknown), mapping only once to reference assembly, flagged in dbSnp as "clinically associated") | 20140222 |
| 1000g2014oct_all       | alternative allele frequency data in 1000 Genomes Project for autosomes | 20141216 |
| esp6500siv2_all        | alternative allele frequency in All subjects in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls. This is lifted over from hg19 by myself. | 20141222 |
| cosmic70              | COSMIC database version 70  | 20140224 |
| clinvar_20150330       | CLINVAR database with Variant Clinical Significance (unknown, untested, non-pathogenic, probable-non-pathogenic, probable-pathogenic, pathogenic, drug-response, histocompatibility, other) and Variant disease name | 20150413 | 
| genomicSuperDups       |  [Duplications of >1000 Bases of Non-RepeatMasked Sequence](http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=440807676_1iJIwAXN34xvNpvISAaGashad4iB&c=chr3&g=genomicSuperDups)  |     20110926    |
| ljb26_all              | whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, PhyloP and SiPhy scores from dbNSFP version 2.6 | 20140925 |
| caddgt10               |     ?                        |    ?    |
| caddindel              |     ?                       |     ?   |
| popfreq_max_20150413   | A database containing all allele frequency from 1000G, ESP6500, ExAC and CG46 | 20150413 |
| mitimpact2             |pathogenicity predictions of human mitochondrial missense variants | 20150520  |

### [selection] [in progress...]
Selection variants of interest:
- FILTER: only keep 'PASS' annotation (remove all variants that did not pass filter criteria)
- INFO: snp138Flagged with rs-number
- INFO: 1000g2014oct_all with AF>0.1
- INFO: ExonicFunc.ensGene = 'synonymous SNV'






