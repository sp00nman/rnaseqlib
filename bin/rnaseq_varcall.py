#!/usr/bin/env python

import argparse
import re
from os import system
import os
import logging



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genetic screen workflow '
                                                 '0.0.1')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False, default="all",
                        choices=["all", "extract", "reorder", "duplicates",
                                 "splitntrim", "bqsr",
                                 "varcall_bamfo", "varcall_samtools",
                                 "varcall_gatk", "filter"],
                        help='Limit job submission to a particular '
                             'Analysis stage.')
    parser.add_argument('--project_name', required=False, type=str,
                        help="name of the project")
    parser.add_argument('--input_file', required=False, type=str,
                        help="name of the sample file")
    parser.add_argument('--sample_dir', required=False, type=str,
                        help="sample directory")
    parser.add_argument('--exec_dir', required=False, type=str,
                        help="exec_dir")
    parser.add_argument('--output_dir', required=False, type=str,
                        help="output_dir")
    parser.add_argument('--ref_genome', required=False, type=str,
                        help="reference genome")
    parser.add_argument('--region', required=False, type=str,
                        help="region eg. 20:30946147-31027122")

    # parse arguments, set defaults
    args = parser.parse_args()
    home_dir = os.getenv("HOME")
    skip_regioncall = False
    if not args.output_dir:
        args.output_dir = os.getcwd()
    if not args.sample_dir:
        args.sample_dir = os.getcwd()
    if not args.exec_dir:
        args.exec_dir = home_dir + "/src"
    if not args.ref_genome:
        args.defuse_ref = home_dir + "/ref_genome"
    if not args.region:
        skip_regioncall = True

    # logging
    # create log file
    logfile_name = args.output_dir + "/" + args.project_name + ".log"
    logging.basicConfig(filename=logfile_name,
                        format='%(levelname)s: %(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)

    # start analysis workflow & logging
    logging.info("RNAseq variant calling (region specific)")

    # start workflow
    if (re.search(r"all|extract", args.stage) and skip_regioncall == False):
        (msg, cmd) = extract(args.input_file, args.sample_dir,
                             args.project_name, args.region, args.output_dir)
        status = run_cmd(msg, cmd)

    if re.search(r"all|reorder", args.stage):
        (msg, cmd) = reorder(args.project_name, args.output_dir,
                             args.ref_genome)
        status = run_cmd(msg, cmd)

    if re.search(r"all|sort_bam", args.stage):
        (msg, cmd) = sort_bam(args.project_name, args.output_dir)
        status = run_cmd(msg, cmd)

    if re.search(r"all|replace_readgroups", args.stage):
        (msg, cmd) = replace_readgroups(args.project_name, args.output_dir)
        status = run_cmd(msg, cmd)

    if re.search(r"all|remove_duplicates", args.stage):
        (msg, cmd) = remove_duplicates(args.project_name, args.output_dir)
        status = run_cmd(msg, cmd)

    if re.search(r"all|index", args.stage):
        (msg, cmd) = index_bam(args.project_name, args.output_dir)
        status = run_cmd(msg, cmd)

    if re.search(r"all|splitntrim", args.stage):
        (msg, cmd) = splitntrim(args.project_name, args.output_dir,
                                args.ref_genome)
        status = run_cmd(msg, cmd)

    if re.search(r"all|bqsr", args.stage):
        (msg, cmd) = bqsr(args.project_name, args.output_dir,
                          args.ref_genome)
        status = run_cmd(msg, cmd)

    if re.search(r"all|varcall_bamfo", args.stage):
        (msg, cmd) = varcall_bamfo(args.project_name, args.output_dir,
                                     args.ref_genome)
        status = run_cmd(msg, cmd)

    if re.search(r"all|varcall_samtools", args.stage):
        (msg, cmd) = varcall_samtools(args.project_name, args.output_dir,
                                     args.ref_genome)
        status = run_cmd(msg, cmd)

    if re.search(r"all|varcall_gatk", args.stage):
        (msg, cmd) = varcall_gatk(args.project_name, args.output_dir,
                                     args.ref_genome)
        status = run_cmd(msg, cmd)

    if re.search(r"all|variant_filtering", args.stage):
        (msg, cmd) = variant_filtering(args.project_name, args.output_dir,
                                       args.ref_genome)
        status = run_cmd(msg, cmd)
