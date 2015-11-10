#!/usr/bin/env python

import argparse
import re
import os
import logging

from rnaseqlib.utils import tools as ts
from rnaseqlib.varcall import annovar_runnables as annovar_rb
from rnaseqlib.varcall import gatk_runnables as gatk
from rnaseqlib.varcall import filter_vcf as fv


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Annotate vcf files for RJ.')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'analysis stage.[process_all,annovar,selection,'
                             'vcf2table]')
    parser.add_argument('--project_name', required=False, type=str,
                        help="name of the project")
    parser.add_argument('--patient_vcf', required=False, type=str,
                        help='List of VCF files/patient to process.')

    parser.add_argument('--sample_dir', required=False, type=str,
                        help="Path to sample directory")
    parser.add_argument('--output_dir', required=False, type=str,
                        help="Path to output directory.")
    # resources
    parser.add_argument('--ref_genome', required=False, type=str,
                        help='Reference genome.')

    # annovar specific options
    parser.add_argument('--annovar', required=False, type=str,
                        help="Annotate variant with annovar.")

    # hardware specific options
    parser.add_argument('--num_cpus', dest='num_cpus', required=False,
                        help='Number of cpus.')
    parser.add_argument('--heap_mem', dest='heap_mem', required=False,
                        help='Maximum heap size provided to Java. [Xmx[num]g]')

    # parse command line arguments
    args = parser.parse_args()

    # set project directory
    project_dir = args.output_dir + "/" + args.project_name

    # create directory structure
    ts.create_output_dir(args.output_dir, args.project_name)

    # load dictionary for file extensions
    data_dir = str(os.path.dirname(os.path.realpath(__file__)).strip('bin')) \
               + "data" + "/"
    file_ext = ts.load_dictionary(data_dir + 'file_extension.txt')
    # load dictionary for stdout messages
    stdout_msg = ts.load_dictionary(data_dir + 'stdout_message.txt')

    #get execution directory
    dn = os.path.dirname(os.path.realpath(__file__))

    # start logging process
    logging.basicConfig(
        filename=args.output_dir + "/" \
                 + args.project_name + "/" \
                 + args.project_name + ".log",
        format='%(levelname)s: %(asctime)s %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG
    )

    # start analysis workflow & logging
    logging.info("Annotate and select vcf files.")

    # start workflow
    if re.search(r"process_all|annovar", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            cmd = annovar_rb.run_annovar(
                vcf_file=path2vcf,
                protocol="ensGene,"
                        + "cytoBand,"
                        + "genomicSuperDups,"
                        + "snp142Mult,"
                        + "snp129,"
                        + "snp142Common,"
                        + "1000g2015feb_all,"
                        + "1000g2015feb_afr,"
                        + "1000g2015feb_amr,"
                        + "1000g2015feb_eas,"
                        + "1000g2015feb_eur,"
                        + "1000g2015feb_sas,"
                        + "esp5400_all,"
                        + "esp6500siv2_all,"
                        + "cosmic70,"
                        + "clinvar_20150330,"
                        + "ljb26_all,"
                        + "caddgt10,"
                        + "caddindel",
                output_file=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id,
                operation="g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f",
                buildversion="hg19",
                annovar_dir=args.annovar
            )

            status = ts.run_cmd(
                message=stdout_msg['annotation'],
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|selection", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            fv.filter_vcf(
                input_file=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id + "."
                        + "hg19_multianno.vcf",
                file_prefix=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id,
                file_suffix=file_ext['inhouse']
            )

    if re.search(r"process_all|vcf2table_proudpv", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            cmd = gatk.var2table_proudpv(
                input_file=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id + "."
                        + file_ext['proud_pv'],
                output_file=project_dir + "/" + args.project_name + "."
                            + uniq_sample_id + "."
                            + file_ext['vcf2table'],
                ref_genome=args.ref_genome,
                heap_mem=args.heap_mem
            )

            status = ts.run_cmd(
                message=stdout_msg['vcf2table'],
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|vcf2table_dna", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            cmd = gatk.var2table_dna(
                input_file=project_dir + "/"
                           + args.project_name + "."
                           + uniq_sample_id + "."
                           + file_ext['proud_pv'],
                output_file=project_dir + "/" + args.project_name + "."
                            + uniq_sample_id + "."
                            + file_ext['vcf2table'],
                ref_genome=args.ref_genome,
                heap_mem=args.heap_mem
            )

            status = ts.run_cmd(
                message=stdout_msg['vcf2table'],
                command=cmd,
                debug=args.debug
            )