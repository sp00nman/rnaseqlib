#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Annotate gene fusions from gene fusion detection tools
- should be independent from tool
- important infromation should be
--> GENE A ; GENE B ; CHR_A; CHR_B; POS_A; POS_B
"""

import sys
import argparse
import re
import logging
import pandas as pd

from rnaseqlib.utils import tools as ts
from rnaseqlib.dragonball import load_fusion_output_files as loadfuse
from rnaseqlib.dragonball import match_fusions as matchf


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Annotate fusions from fusion tools for RNA-seq.')
    parser.add_argument(
        '--project_name', required=False, type=str,
        help="Name of the project directory.")
    parser.add_argument(
        '--input_file', required=False, type=str,
        help="Input file with the format (tab separated): "
             "tool[defuse,tophatfusion,soapfuse] path/to/file")
    parser.add_argument(
        '--annotation_file', required=False, type=str,
        help="Annotation file with the format (tab separated): "
             "database[see github] /path/to/annotation_file [pairs(A_B)/individual(A,B) [annotate/filter]")
    parser.add_argument(
        '--output_dir', required=False, type=str,
        hep="Path to output directory.")

    # hardware specific options
    parser.add_argument(
        '--num_cpus', dest='num_cpus', required=False, help='Number of cpus.')

    # parse command line arguments
    args = parser.parse_args()

    # set project directory
    project_dir = args.output_dir + "/" + args.project_name
    ts.create_output_dir(args.output_dir, args.project_name)

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

    fusion_files = ts.load_tab_delimited(args.input_file)
    annotation = pd.read_csv(args.annotation_file, sep="\t")

    # for each annotation run all samples !!!!
    annotation = pd.read_csv(
        args.annotation_file,
        sep="\t")

    samples = pd.read_csv(
        args.input_file,
        sep="\t",
        names=["tool", "path"])

    for sample in samples.iterrows():

        sample_sr = pd.Series(sample)
        tool = sample_sr.get(1)['tool']
        path_sample = sample_sr.get(1)['path']

        fusions = loadfuse.load_fusion(path_sample, tool)

        for annotation in annotation.iterrows():

            # for each annotation...
            annotation_sr = pd.Series(annotation)
            label = annotation_sr.get(1)['label']
            gene_position = annotation_sr.get(1)['gene_position']
            pair_single_distance = annotation_sr.get(1)['pair_single_distance']
            source = annotation_sr.get(1)['source']
            filter_annotate = annotation_sr.get(1)['filter_annotate']
            file_location = annotation_sr.get(1)['file_location']

            if pair_single_distance is "pair":
                fusions[annotation_sr] = fusions.apply(
                    lambda row: matchf.filter_gene_pair(
                        row['gene_A'],
                        row['gene_B']
                    ),
                    axis=1
                )



