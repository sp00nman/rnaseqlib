### rnaseq_varcall.py -h

```bash
usage: rnaseq_varcall.py [-h] [--debug DEBUG]
                         [--stage {all,extract,reorder,duplicates,splitntrim,realignment,recalibration,varcall,filter}]
                         [--project_name PROJECT_NAME]
                         [--input_file INPUT_FILE] [--sample_dir SAMPLE_DIR]
                         [--exec_dir EXEC_DIR] [--output_dir OUTPUT_DIR]
                         [--ref_genome REF_GENOME] [--region REGION]

Genetic screen workflow 0.0.1

optional arguments:
  -h, --help            show this help message and exit
  --debug DEBUG         Debug level
  --stage {all,extract,reorder,duplicates,splitntrim,realignment,recalibration,varcall,filter}
                        Limit job submission to a particular Analysis stage.
  --project_name PROJECT_NAME
                        name of the project
  --input_file INPUT_FILE
                        name of the sample file
  --sample_dir SAMPLE_DIR
                        sample directory
  --exec_dir EXEC_DIR   exec_dir
  --output_dir OUTPUT_DIR
                        output_dir
  --ref_genome REF_GENOME
                        reference genome
  --region REGION       region eg. 20:30946147-31027122
´´´
