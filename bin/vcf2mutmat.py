#!/usr/bin/env python

import argparse
import re
import os
import logging
from rnaseqlib.utils import tools as ts
from rnaseqlib.varcall import gatk_runnables as gatk
from rnaseqlib.varcall import extract_genes_of_interest as goi


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert vcf files to patient '
                                                 'mutation matrix')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'analysis stage.[all,vcf2table,targetlist,'
                             'plot]')
    parser.add_argument('--output_dir', required=False, type=str,
                        help="Project directory output.")
    parser.add_argument('--project_name', required=False, type=str,
                        help="name of the project")
    parser.add_argument('--patient_vcf', required=False, type=str,
                        help="List of VCF files/patient to process.")
    parser.add_argument('--target_list', required=False, type=str,
                        help="List of genes of interest.")
    parser.add_argument('--ref_genome', required=False, type=str,
                        help="reference genome")
    parser.add_argument('--id_conversion', required=False, type=str,
                        help="Conversion table.")
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
    logging.info("Convert vcf files to patient mutation matrix.")

    # start workflow
    if re.search(r"all|vcf2table", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            cmd = gatk.variants2table(
                input_file=path2vcf,
                output_file=project_dir + "_"
                            + uniq_sample_id
                            + ".vcf2table",
                ref_genome=args.ref_genome,
                heap_mem=args.heap_mem
            )

            status = ts.run_cmd(
                message=stdout_msg['var2table'],
                command=cmd,
                debug=args.debug
            )

    if re.search(r"all|targetlist", args.stage):

        goi.extract_genes(
            project_dir=project_dir,
            output_file=project_dir + "_"
                         "patient_mut.table",
            target_list=args.target_list,
            vcf2table=args.patient_vcf,
            conversion_table=args.id_conversion
        )

    if re.search(r"all|plot", args.stage):
        pass
