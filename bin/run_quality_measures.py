#!/usr/bin/env python

import argparse
import re
import os
import logging

from rnaseqlib.quality import rseqc_runnables as rseqc
from rnaseqlib.quality import fastqc_runnables as fastqc
from rnaseqlib.utils import tools as ts

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Run quality measures for RNA-seq on BAM files')
    parser.add_argument(
        '--debug', dest='debug', required=False, type=int, help='Debug level')
    parser.add_argument(
        '--stage', dest='stage', required=False, help='Limit job submission '
        'to a particular analysis stage.[bam_stat,genebody_coverage,'
        'inner_distance, junction_annotation, junction_saturation, '
        'read_duplication, read_gc, fastqc, process_all]')
    parser.add_argument('--project_name', required=False, type=str,
                        help="name of the project")
    parser.add_argument('--patient_bam', required=False, type=str,
                        help='List of BAM files/patient to process.')
    parser.add_argument('--output_dir', required=False, type=str,
                        help="Path to output directory.")
    # resources
    parser.add_argument('--ref_genome', required=False, type=str,
                        help='Reference genome in BED format.')

    # hardware specific options
    parser.add_argument('--num_cpus', dest='num_cpus', required=False,
                        help='Number of cpus.')
    parser.add_argument('--heap_mem', dest='heap_mem', required=False,
                        help='Maximum heap size provided to Java. [Xmx[num]g]')

    # parse command line arguments
    args = parser.parse_args()

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
    logging.info("Run quality measures on RNA-seq BAM files.")

    # load bamfiles
    bam_files = ts.load_tab_delimited(args.patient_bam)

    if re.search(r"process_all|bam_stat", args.stage):

        for bam in bam_files:

            uniq_sample_id = bam[0]
            path2bam = bam[2]

            cmd = rseqc.bam_stat(
                path2bam,
                "rseqc_bam_stat" + "." + uniq_sample_id + ".txt"
            )

            status = ts.run_cmd(
                message="bam_stat.py started",
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|genebody_coverage", args.stage):

        for bam in bam_files:

            uniq_sample_id = bam[0]
            path2bam = bam[2]

            cmd = rseqc.genebody_coverage(
                path2bam,
                args.reference ,
                "rseqc_genebody_coverage" + "." + uniq_sample_id + ".txt"
            )

            status = ts.run_cmd(
                message="genebody_coverage.py started",
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|inner_distance", args.stage):

        for bam in bam_files:

            uniq_sample_id = bam[0]
            path2bam = bam[2]

            cmd = rseqc.inner(
                path2bam,
                args.reference,
                "reseqc_inner_distance" + "." + uniq_sample_id + ".txt"
            )

            status = ts.run_cmd(
                message="inner_distance.py started",
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|junction_annotation", args.stage):

        for bam in bam_files:

            uniq_sample_id = bam[0]
            path2bam = bam[2]

            cmd = rseqc.genebody_coverage(
                path2bam,
                args.reference,
                "reseqc_inner_distance" + "." + uniq_sample_id + ".txt"
            )

            status = ts.run_cmd(
                message="inner_distance.py started",
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|junction_saturation", args.stage):

        for bam in bam_files:

            uniq_sample_id = bam[0]
            path2bam = bam[2]

            cmd = rseqc.junction_saturation(
                path2bam,
                args.reference,
                "reseqc_junction_saturation" + "." + uniq_sample_id + ".txt"
            )

            status = ts.run_cmd(
                message="junction_saturation.py started",
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|read_duplication", args.stage):

        for bam in bam_files:

            uniq_sample_id = bam[0]
            path2bam = bam[2]

            cmd = rseqc.read_duplication(
                path2bam,
                "reseqc_read_duplication" + "." + uniq_sample_id + ".txt"
            )

            status = ts.run_cmd(
                message="inner_distance.py started",
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|read_gc", args.stage):

        for bam in bam_files:

            uniq_sample_id = bam[0]
            path2bam = bam[2]

            cmd = rseqc.read_gc(
                path2bam,
                "reseqc_read_gc" + "." + uniq_sample_id + ".txt"
            )

            status = ts.run_cmd(
                message="read_gc.py started",
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|fastqc", args.stage):

        for bam in bam_files:

            uniq_sample_id = bam[0]
            path2bam = bam[2]

            cmd = fastqc.fastqc(
                path2bam,
                args.output_dir,
            )

            status = ts.run_cmd(
                message="fastqc started",
                command=cmd,
                debug=args.debug
            )
