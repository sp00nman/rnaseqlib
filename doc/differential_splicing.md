
Differential splicing analysis for RNA-seq
==========================================
(in preparation)

Differential splicing analysis was not implemented as a workflow, instead it requires several independent steps.

#### STEP 1
Mapping reads with [STAR](https://github.com/alexdobin/STAR). 
This produces output files named *.SJ.out.tab which will be used as input for STEP 2

#### STEP 2


```
#!/bin/bash
#SBATCH --job-name=[job_name]
#SBATCH --cpus-per-task=[32]
#SBATCH --distribution=[block]
#SBATCH --mem=[256000]
#SBATCH --time=[14-00:00:00]
#SBATCH --partition=[longq]
#SBATCH --account=[user_name]
#SBATCH -o [output_log_file]
#SBATCH -e [error_log_file]

Rscript --vanilla \
bin/differential_splicing.R \
[meta_data]
[exp_name]
[design_formulas]
[working_dir]
[project_name]
[num_cores]

/scratch/lab_kralovics/fschischlik/prj_psi/meta_data_PMFControlonly.tsv \
simple_model \
/scratch/lab_kralovics/fschischlik/prj_psi/design_formulas.txt \
/scratch/lab_kralovics/fschischlik/prj_psi/exp12_numpat70_jxnfilt10 \
exp12_numpat70_jxnfilt10 \
32
```
