#!/usr/bin/env python
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
from rnaseqlib.utils import convert_gene_ids as cv


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
        help="Path to output directory.")
    parser.add_argument(
        '--id_conversion', required=False, type=str,
        help='Conversion table. Ensembl ids to genesymbols')

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
    logging.info("Annotate fusions from fusion tools for RNA-seq.")

    fusion_files = ts.load_tab_delimited(args.input_file)
    annotation = pd.read_csv(args.annotation_file, sep="\t")

    samples = pd.read_csv(
        args.input_file,
        sep="\t",
        names=["uniq_sample_id", "tool", "path"])

    # load resources
    ens2gs_conversion_table = cv.read_ensgene_genesymb(
        args.id_conversion
    )

    for sample in samples.iterrows():

        sample_sr = pd.Series(sample)
        sample = sample_sr.get(1)['uniq_sample_id']
        tool = sample_sr.get(1)['tool']
        path_sample = sample_sr.get(1)['path']

        print sample
        print tool
        print path_sample

        fusions = loadfuse.load_fusion(path_sample, tool)

        # need to load again and again otherwise I get an
        # error message....duuuhhh
        annotation = pd.read_csv(
            args.annotation_file,
            sep="\t")

        for annotation in annotation.iterrows():
            # for each annotation...
            annotation_sr = pd.Series(annotation)
            label = annotation_sr.get(1)['label']
            gene_position = annotation_sr.get(1)['gene_position']
            pair_single_distance = annotation_sr.get(1)['pair_single_distance']
            source = annotation_sr.get(1)['source']
            filter_annotate = annotation_sr.get(1)['filter_annotate']
            file_location = annotation_sr.get(1)['file_location']

            #prints
            print label
            print gene_position
            print pair_single_distance
            print source
            print file_location
            print filter_annotate
            label_name = label + "_" + source
            # read in annotation file
            annotation_dta = ts.read_annotation_file(file_location)

            if pair_single_distance == "pair":

                fusions[label_name] = fusions.apply(
                    lambda row: matchf.annotate_gene_pair(
                        row['gene_A'],
                        row['gene_B'],
                        ens2gs_conversion_table,
                        annotation_dta,
                        tool
                    ),
                    axis=1
                )

            if pair_single_distance == "single":

                fusions[label_name] = fusions.apply(
                    lambda row: matchf.annotate_single(
                        row['gene_A'],
                        row['gene_B'],
                        ens2gs_conversion_table,
                        annotation_dta,
                        tool
                    ),
                    axis=1
                )

        output_file = project_dir + "/" \
                      + args.project_name + "." \
                      + sample + "_" \
                      + "annotated.txt"

        out_handle = open(output_file, 'w')

        fusions.to_csv(
            out_handle,
            sep="\t",
            index=0,
            header=True
        )

        out_handle.close()

