#!/usr/bin/env python

import argparse
import re
import os
import logging
from varcall import runnables as rb
from utils import tools as ts


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genetic screen workflow '
                                                 '0.0.1')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False, default="all",
                        choices=["all", "alignment", "region", "duplicates",
                                 "splitntrim", "bqsr", "varcall_bamfo",
                                 "varcall_samtools", "varcall_gatk", "filter"],
                        help='Limit job submission to a particular '
                             'analysis stage.')
    parser.add_argument('--project_name', required=False, type=str,
                        help="name of the project")
    parser.add_argument('--read1', required=False, type=str,
                        help="For paired alignment, forward read.")
    parser.add_argument('--read2', required=False, type=str,
                        help="For paired alignment, reverse read.")
    parser.add_argument('--sample_dir', required=False, type=str,
                        help="sample directory")
    parser.add_argument('--exec_dir', required=False, type=str,
                        help="exec_dir")
    parser.add_argument('--output_dir', required=False, type=str,
                        help="output_dir")
    parser.add_argument('--ref_genome', required=False, type=str,
                        help="reference genome")
    parser.add_argument('--region', required=False, type=str, default=False,
                        help="region eg. 20:30946147-31027122")
    parser.add_argument('--num_cpus', dest='num_cpus', required=False,
                        help='Number of cpus.')

    # defaults
    args = parser.parse_args()
    home_dir = os.getenv("HOME")

    if not args.output_dir:
        args.output_dir = os.getcwd()
    if not args.sample_dir:
        args.sample_dir = os.getcwd()
    if not args.exec_dir:
        args.exec_dir = home_dir + "/src"
    if not args.ref_genome:
        args.defuse_ref = home_dir + "/ref_genome"
    if not args.num_cpus:
        args.num_cpus = "1"

    # create log file
    logfile_name = args.output_dir + "/" + args.project_name + "/" \
                   + args.project_name + ".log"
    logging.basicConfig(filename=logfile_name,
                        format='%(levelname)s: %(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)

    # start analysis workflow & logging
    logging.info("RNAseq variant calling (region specific)")

    # load dictionary for file extensions
    data_dir = str(os.path.dirname(os.path.realpath(__file__)).strip('bin')) \
               + "data" + "/"
    file_ext = ts.load_dictionary(data_dir + 'file_extension.txt')
    # load dictionary for stdout messages
    stdout_msg = ts.load_dictionary(data_dir + 'stdout_message.txt')

    # start workflow
    if re.search(r"all|alignment", args.stage):
        cmd = rb.rnaseq_align(
            genome_path=args.genome,
            read1=args.read1,
            read2=args.read2,
            num_cpus=args.num_cpus
        )

        status = ts.run_cmd(
            message=stdout_msg['alignment'],
            command=cmd,
            debug=args.debug
        )

        cmd = rb.star_index(
            genome_path=args.ref_genome,
            genome="hg19.fa",
            firstroundalignment="dummy",
            sjdbOverhang="75",
            num_cpus=args.num_cpus
        )

        status =ts.run_cmd(
            message=stdout_msg['alignment_index'],
            command=cmd,
            debug=args.debug
        )

        cmd = rb.rnaseq_align(
            genome_path=different_path,
            read1=args.read1,
            read2=args.read2,
            num_cpus=args.num_cpus
        )

    if re.search(r"all|extract", args.stage) and args.skip_regioncall is False:
        cmd = rb.extract(
            input_file=args.input_file,
            region=args.sample_dir,
            output_file=args.output_dir)

        status = ts.run_cmd(
            message=stdout_msg['extract'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|reorder", args.stage):
        cmd = rb.reorder_sam(
            inbamfile=args.project_name,
            outbamfile=args.output_dir,
            genome_path=args.ref_genome)

        status = ts.run_cmd(
            message=stdout_msg['reorder'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|sort_bam", args.stage):
        cmd = rb.sort_bam(
            inbamfile=args.project_name,
            outbamfile=args.output_dir,
            sort_order="coordinate")

        status = ts.run_cmd(
            message=stdout_msg['sort_bam'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|replace_readgroups", args.stage):
        cmd = rb.replace_readgroups(
            input_file=args.project_name,
            output_file=args.output_dir,
            project_name=args.project_name)

        status = ts.run_cmd(
            message=stdout_msg['replace_rg'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|remove_duplicates", args.stage):
        cmd = rb.remove_duplicates(
            inbamfile=args.project_name,
            outbamfile=args.output_dir,
            metrics_file=".duplicate_metrics.txt")

        status = ts.run_cmd(
            message=stdout_msg['duplicates'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|index", args.stage):
        cmd = rb.index_bam(
            inbamfile=args.project_name)

        status = ts.run_cmd(
            message=stdout_msg['index'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|splitntrim", args.stage):
        cmd = rb.splitntrim(
            input_file=args.project_name,
            output_file=args.output_dir,
            ref_genome=args.ref_genome)

        status = ts.run_cmd(
            message=stdout_msg['splitntrim'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|bqsr", args.stage):
        cmd = rb.bqsr(
            input_file=args.project_name,
            output_file=args.output_dir,
            ref_genome=args.ref_genome)

        status = ts.run_cmd(
            message=stdout_msg['bqsr'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|varcall_bamfo", args.stage):
        cmd = rb.varcall_bamfo(
            input_file=args.project_name,
            output_file_gatk=args.output_dir,
            ref_genome=args.ref_genome)

        status = ts.run_cmd(
            message=stdout_msg['bamfo'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|varcall_samtools", args.stage):
        cmd = rb.varcall_samtools(
            input_file=args.project_name,
            output_file_gatk=args.output_dir,
            ref_genome=args.ref_genome)

        status = ts.run_cmd(
            message=stdout_msg['samtools'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|varcall_gatk", args.stage):
        cmd = rb.varcall_gatk(
            input_file=args.project_name,
            output_file_gatk=args.output_dir,
            ref_genome=args.ref_genome)

        status = ts.run_cmd(
            message=stdout_msg['gatk'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|variant_filtering", args.stage):
        cmd = rb.variant_filtering(
            input_file=args.project_name,
            output_file=args.output_dir,
            ref_genome=args.ref_genome)

        status = ts.run_cmd(
            message=stdout_msg['filtering'],
            command=cmd,
            debug=args.debug
        )