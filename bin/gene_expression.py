#!/usr/bin/env python

import argparse
import re
import os
import logging

from rnaseqlib.utils import tools as ts
from rnaseqlib.varcall import star_runnables as star_rb
from rnaseqlib.expression import htseq_runnables as htseq
from rnaseqlib.varcall import samtools_runnables as samtools

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='RNA-seq gene-expression')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'analysis stage.[process_all,sort,htseq]')
    parser.add_argument('--project_name', required=False, type=str,
                        help="name of the project")
    parser.add_argument('--sample_list', required=False, type=str,
                        help='List of samples to process.')

    parser.add_argument('--output_dir', required=False, type=str,
                        help="Path to output directory.")

    parser.add_argument('--gtf', required=False, type=str,
                        help="gtf file annotation")
    parser.add_argument('--stranded', required=False, type=str,
                        help="Is the library stranded? [yes,no]")

    # parse command line arguments
    args = parser.parse_args()

    # set project directory
    project_dir = args.output_dir + "/" + args.project_name

    # create directory structure
    ts.create_output_dir(args.output_dir, args.project_name)
    sub_dir = args.output_dir + "/" + args.project_name

    # create log file
    logfile_name = args.output_dir + "/" + args.project_name + "/" \
                   + args.project_name + ".log"
    logging.basicConfig(filename=logfile_name,
                        format='%(levelname)s: %(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)

    # start analysis workflow & logging
    logging.info("Initiated RNA-seq expression analysis workflow")

    # load dictionary for file extensions
    data_dir = str(os.path.dirname(os.path.realpath(__file__)).strip('bin')) \
               + "data" + "/"
    file_ext = ts.load_dictionary(data_dir + 'file_extension.txt')
    # load dictionary for stdout messages
    stdout_msg = ts.load_dictionary(data_dir + 'stdout_message.txt')

    # start workflow
    if re.search(r"process_all|sort", args.stage):

        bamfiles = ts.load_tab_delimited(args.sample_list)

        for bamfile in bamfiles:

            uniq_sample_id = bamfile[0]
            patient_id = bamfile[1]
            path2bam = bamfile[2]

            cmd = samtools.sort_bam(
                inbamfile=path2bam,
                output_file=project_dir + "/"
                        + uniq_sample_id + "_"
                        + file_ext['simple_sort']
            )

            status = ts.run_cmd(
                message=stdout_msg['sort'],
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|count", args.stage):

        bamfiles = ts.load_tab_delimited(args.sample_list)

        for bamfile in bamfiles:

            uniq_sample_id = bamfile[0]
            patient_id = bamfile[1]
            path2bam = bamfile[2]

            cmd = htseq.htseq_count(
                uniq_sample_id=uniq_sample_id,
                output_dir=project_dir,
                bamfile=project_dir + "/"
                        + uniq_sample_id + "_"
                        + file_ext['simple_sort'],
                outfile=project_dir + "/"
                        + uniq_sample_id + "_"
                        + file_ext['counts'],
                stranded=args.stranded,
                order="name",
                minaqual="10",
                mode="union",
                gtf=args.gtf
            )

            status = ts.run_cmd(
                message=stdout_msg['htseq'],
                command=cmd,
                debug=args.debug
            )