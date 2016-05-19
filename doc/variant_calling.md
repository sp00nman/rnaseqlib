### Workflow for calling variants optimized for RNA sequencing data.

The following steps are well described [here] (http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq)
[alignment,star2pass,replace_rg,duplicates,index,splitntrim,indel,bqsr,gatk, gatk_flag]
Part of the workflow are in-house implementations:
[hrun_flag, nind_flag,annovar,selection]


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
                        bamfo,samtools,gatk,gatk_flag,hrun_flag,nind_flag,annovar,selection]
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

## Variant Annotation & Filtering Steps

### [gatk_flag]
GATK ```-T VariantFiltration``` is used to flag the following variants:
- [1] at least 3 SNPs that are within a window of 35 bases ```-window 35 -cluster 3```
 + ```##FILTER=<ID=SnpCluster,Description="SNPs found in clusters">```
- [2] fisher strand value FS>30
 + ```##FILTER=<ID=FS,Description="FS > 30.0">```
- [3] Quality by Depth QD<2
 + ```##FILTER=<ID=QD,Description="QD < 2.0">```


### [hrun_flag]
- [4] flag variants within homopolymer runs >=5
 + ```[##FILTER=<ID=HRun,Description="HRun >= 5">]```
- [5] flag variants within 1bp away from a homopolymer runs >=5; nHrun (near homopolymer run)
 + ```[##FILTER=<ID=nHRun,Description="nHRun >= 5">]```

### [nind_flag] 
- [6] flag SNVs located within 5 bases from an indel (Reumers et al., Nature Biotech., 2012)
 + ```[##FILTER=<ID=nIndel,Description="nIndel <= 5">]```

### [annovar]
Databases used for annotation [ANNOVAR] (http://annovar.openbioinformatics.org/en/latest/user-guide/download/) and
[inhouse] filtering steps:

| NAME                   | DESCRIPTION                 | DATE   |
| :--------------------- |:----------------------------|:-------|
| cytoBand               | cytoBand annotation             |   ?     |
| genomicSuperDups  |  [Duplications of >1000 Bases of Non-RepeatMasked Sequence](http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=440807676_1iJIwAXN34xvNpvISAaGashad4iB&c=chr3&g=genomicSuperDups)  |     20110926    |
| snp142Mult              | ?	| ? |
| snp129          | dbSNP129 with dbSNP with ANNOVAR index files | 20120809 |
| snp142Common    | dbSNP142 only with common SNPs, MAF>0.1; downloaded from ucsc without dbSNP with ANNOVAR index files; index is produced on the fly| 20150822 |
| 1000g2015feb_all       | alternative allele frequency data in 1000 Genomes Project for autosomes. | 201502?? |
| 1000g2015feb_afr       | same as 1000g2015feb_all for AFR (African) [691 total] | 201502?? |
| 1000g2015feb_amr       | same as 1000g2015feb_all for AMR (Admixed American) [355 total] | 201502?? |
| 1000g2015feb_eas       | same as 1000g2015feb_all for EAS (East Asian) [523 total]  | 201502?? |
| 1000g2015feb_eur       | same as 1000g2015feb_all for EUR (European) [514 total] | 201502?? |
| 1000g2015feb_sas       | same as 1000g2015feb_all for SAS (South Asian) [494 total] | 201502?? |
| esp5400_all        | [Exome Variant Server] alternative allele frequency in All subjects in the NHLBI-ESP project with 5400 exomes. This is lifted over from hg19. (Includes sequencing errors.) | ?????? |
| esp6500siv2_all        | [Exome Variant Server] alternative allele frequency in All subjects in the NHLBI-ESP project with 6500 exomes (2203 African-Americans and 4300 European-Americans unrelated individuals), including the indel calls and the chrY calls. This is lifted over from hg19. | 20141222 |
| cosmic70              | COSMIC database version 70  | 20140224 |
| clinvar_20150330       | CLINVAR database with Variant Clinical Significance (unknown, untested, non-pathogenic, probable-non-pathogenic, probable-pathogenic, pathogenic, drug-response, histocompatibility, other) and Variant disease name | 20150413 | 
| ljb26_all              | whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, PhyloP and SiPhy scores from dbNSFP version 2.6 | 20140925 |
| caddgt10               |     ?                        |    ??????    |
| caddindel              |     ?                       |     ??????   |

### [selection] [in progress...]
Selection variants of interest (somatic & cancer relevant mutations):
- FILTER: only keep 'PASS' annotation (remove all variants that did not pass filter criteria)
- INFO: snp142Common with rs-number
- INFO: 1000g2014oct_all with MAF>0.01
- INFO: esp5400siv2_all with MAF>0.01
- INFO: ExonicFunc.ensGene = 'synonymous SNV'






