### rnaseq_varcall.py -h

```bash
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
                        ll,alignment,extract,duplicates,splitntrim,bqsr,bamfo,
                        samtools,gatk,filter,annotation]
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
´´´
