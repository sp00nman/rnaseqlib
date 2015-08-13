### Workflow for calling variants optimized for RNA sequencing data.

Workflow is well described [here] (http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq)

## Software requirements
 
+ [STAR] (https://github.com/alexdobin/STAR)
+ [Samtools] (http://samtools.sourceforge.net/)
+ [GATK] (https://www.broadinstitute.org/gatk/)
+ [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)

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
                        bamfo,samtools,gatk,filter,annotation]
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
If ```--star2pass``` is set then a new index is created using splice junction information contained in the file SJ.out.tab from the first pass.

### [replace_rg]
Adds readgroup information and sorts by coordinate as described [here] (http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq)

### [duplicates]
Remove duplicate reads as described [here] (http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq)

### [splitntrim]
GATK tool ```-T SplitNCigarReads``` splits reads into exon segments and hard-clip sequences overhanging into the
intronic regions. Additionally mapping qualities are reassigned.

### [bqsr]
Not implemented yet.

### [gatk]
GATK tool ```-T HaplotypeCaller``` caller with options.
```bash
-dontUseSoftClippedBases
-recoverDanglingHeads
-stand_call_conf 20.0

```
Minimum phred-scaled confidence is lowered to 20 according to GATK best practices recommendation.

### [gatk-filter]
GATK ```-T VariantFiltration``` is used to filter for:
- [1] at least 3 SNPs that are within a window of 35 bases ```-window 35 -cluster 3```
- [2] fisher strand value FS>30
- [3] Quality by Depth QD<2

### [inhouse-filter] [not implemented yet]
- [4] filter for variants within homopolymer runs >=5 ```[##FILTER=<ID=HRun,Description="HRun >= 5">]```
- [5] remove sites from repetitive regions (according to RepeatMasker annotation)
- [6] filter for common variants in snp138NonFlagged

#### suggested (Reumers et al., Nature Biotech., 2012):
- [7] filter for SNVs located within 5 bases from and indel

### [annotation]
*.vcf files are converted to ANNOVAR file format and filtered for common variant that appear in snp138NonFlagged




