```
usage: vcf2mutmat.py [-h] [--debug DEBUG] [--stage STAGE]
                     [--output_dir OUTPUT_DIR] [--project_name PROJECT_NAME]
                     [--patient_vcf PATIENT_VCF] [--target_list TARGET_LIST]
                     [--ref_genome REF_GENOME] [--id_conversion ID_CONVERSION]
                     [--heap_mem HEAP_MEM]

Convert vcf files to patient mutation matrix

optional arguments:
  -h, --help            show this help message and exit
  --debug DEBUG         Debug level
  --stage STAGE         Limit job submission to a particular analysis
                        stage.[all,vcf2table,targetlist,plot]
  --output_dir OUTPUT_DIR
                        Project directory output.
  --project_name PROJECT_NAME
                        name of the project
  --patient_vcf PATIENT_VCF
                        List of VCF files/patient to process.
  --target_list TARGET_LIST
                        List of genes of interest.
  --ref_genome REF_GENOME
                        reference genome
  --id_conversion ID_CONVERSION
                        Conversion table.
  --heap_mem HEAP_MEM   Maximum heap size provided to Java. [Xmx[num]g]

```
