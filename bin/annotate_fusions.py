#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Annotate gene fusions from gene fusion detection tools
- should be independent from tool
- important infromation should be
--> GENE A ; GENE B ; CHR_A; CHR_B; POS_A; POS_B
"""

import argparse
import logging
import pandas as pd
import gffutils
import pyBigWig
import pybedtools

from rnaseqlib.utils import tools as ts
from rnaseqlib.utils import convert_gene_ids as cv
from rnaseqlib.dragonball import load_fusion_output_files as loadfuse
from rnaseqlib.dragonball import match_fusions as matchf
from rnaseqlib.dragonball import gene_gene_distance as ggdis
from rnaseqlib.dragonball import get_biotype as gbio
from rnaseqlib.dragonball import get_mapability as getmap
from rnaseqlib.dragonball import filter_fusions as ffus


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
    # resources
    parser.add_argument(
        '--gtf_dbfile', required=False, type=str,
        help='gtf file (transformed to database with gffutils)')
    parser.add_argument(
        '--id_conversion', required=False, type=str,
        help='Conversion table. Ensembl ids to genesymbols')

    # hard cut off filters
    parser.add_argument(
        '--min_num_split_reads', required=False, type=str,
        help="Minimum amount of split reads.")
    parser.add_argument(
        '--min_num_span_reads', required=False, type=str,
        help="Minimum amount of spanning reads."
    )

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

    ## load resources ##

    # load fusion sample file
    #fusion_files = ts.load_tab_delimited(args.input_file)
    samples = pd.read_csv(
        args.input_file,
        sep="\t",
        names=["uniq_sample_id", "tool", "path"])

    # load annotations
    annotation = pd.read_csv(
        args.annotation_file, sep="\t")

    # load id conversion table
    ens2gs_conversion_table = cv.read_ensgene_genesymb(
        args.id_conversion)

    # load gtf database file
    # this is done at this step, in order to avoid reading it again and
    # again for each sample....
    gtfdb = gffutils.FeatureDB(
        args.gtf_dbfile, keep_order=True)

    gene_version = {}

    for gene in gtfdb.features_of_type('gene'):
        # gencode version 7 uses ensembl ids with version
        # number for each gene/transcript (eg. ENSG00000096968.8)
        # need to create a hash table that stores the "raw"
        # ensembl id together with the gene/transcript version
        # number; this will later be used to query for that gene
        # within the database...sooo stupid...
        primary_key = gene.id
        (ensembl, version_num) = primary_key.split(".")
        gene_version[ensembl] = version_num

    ## start iterations for each sample ##
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
            mode = annotation_sr.get(1)['mode']
            source = annotation_sr.get(1)['source']
            filter_annotate = annotation_sr.get(1)['filter_annotate']
            file_location = annotation_sr.get(1)['file_location']

            label_name = label + "_" + source

            print label_name

            if mode == "pair":

                # read annotation file
                annotation_dta = ts.read_annotation_file(file_location)

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

            if mode == "single":

                # read annotation file
                annotation_dta = ts.read_annotation_file(file_location)

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

            if mode == "distance":

                fusions[label_name] = fusions.apply(
                    lambda row: ggdis.min_dis_gene(
                        row['gene_A'],
                        row['gene_B'],
                        ens2gs_conversion_table,
                        tool,
                        gtfdb,
                        gene_version
                    ),
                    axis=1
                )

            if mode == "biotype":

                if label == "5_prime_biotype":
                    gene = 'gene_A'
                else:
                    gene = 'gene_B'

                fusions[label_name] = fusions.apply(
                    lambda row: gbio.get_biotype(
                        row[gene],
                        gtfdb,
                        tool,
                        ens2gs_conversion_table,
                        gene_version
                    ),
                    axis=1
                )
            if mode == "mapability_sequence_pos":

                # read in bigwig file
                bw = pyBigWig.open(file_location)

                if label == "5_prime_mapability":
                    chrom = 'chr_A'
                    breakpoint = 'genomicbreakpoint_A'
                else:
                    chrom = 'chr_B'
                    breakpoint = 'genomicbreakpoint_B'

                fusions[label_name] = fusions.apply(
                    lambda row: getmap.get_mean_values(
                        bw,
                        row[chrom],
                        row[breakpoint]
                    ),
                    axis=1
                )

            if mode == "blacklisted_regions":

                # read in bedfile
                bedfile = pybedtools.BedTool(file_location)

                if label == "5_prime_exclude":
                    chrom = 'chr_A'
                    breakpoint = 'genomicbreakpoint_A'
                else:
                    chrom = 'chr_B'
                    breakpoint = 'genomicbreakpoint_B'

                fusions[label_name] = fusions.apply(
                    lambda row: getmap.intersect_genomic_breakpoint(
                        bedfile,
                        row[chrom],
                        row[breakpoint]
                    ),
                    axis=1
                )
        # add column FILTER:
        fusions['FILTER'] = fusions.apply(
            lambda row: ffus.match_citeria(row), axis=1)

        output_file = project_dir + "/" \
                      + args.project_name + "." \
                      + sample + "_" + tool + "_"\
                      + "annotated.txt"

        out_handle = open(output_file, 'w')

        fusions.to_csv(
            out_handle,
            sep="\t",
            index=0,
            header=True
        )

        out_handle.close()

