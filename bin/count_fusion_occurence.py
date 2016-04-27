#!/usr/bin/env python

import argparse
import re
import os
import logging
from rnaseqlib.utils import tools as ts
from rnaseqlib.varcall import gatk_runnables as gatk
from rnaseqlib.utils import convert_gene_ids as cv
import pandas as pd

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Count number of occurrences of fusions in samples.')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'analysis stage.[all]')
    parser.add_argument('--output_dir', required=False, type=str,
                        help='Project directory output.')
    parser.add_argument('--project_name', required=False, type=str,
                        help='Name of the project.')
    parser.add_argument('--input_list', required=False, type=str,
                        help='List of files to include')

    # resources
    parser.add_argument('--id_conversion', required=False, type=str,
                        help='Conversion table. Ensembl ids to genesymbols')

    # parse command line arguments
    args = parser.parse_args()

    # set project directory
    project_dir = args.output_dir + "/" + args.project_name

    # create directory structure
    ts.create_output_dir(args.output_dir, args.project_name)

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
    logging.info("Count number of occurrences of fusions in samples.")

    if re.search(r"all|append_fusions", args.stage):

        output_file = project_dir + "/" + args.project_name + "_" \
                      + "fusion_table.txt"
        out_handle = open(output_file, 'a')

        fusion_file = ts.load_tab_delimited(args.input_list)

        for file in fusion_file:

            uniq_sample_id = file[0]
            patient_id = file[1]
            file_path = file[2]

            fusion_table = pd.read_csv(
                file_path,
                sep="\t",
                names=["RANKING",
                        "DETECTION_SOFTWARE",
                        "UNIQ_SAMPLE_ID",
                        "gene_A",
                        "gene_B",
                        "chr_A",
                        "chr_B",
                        "genomicbreakpoint_A",
                        "genomicbreakpoint_B",
                        "splitr_count",
                        "span_count",
                        "score",
                        "FILTER"],
                header=True
            )

            # add patient information
            fusion_table['UNIQ_SAMPLE_ID_ANNO'] = uniq_sample_id
            fusion_table['PATIENT_ID'] = patient_id

            # convert genesymbols to ensids
            key_file = cv.read_ensgene_genesymb(
                args.id_conversion,
                mode="genesymbol2ensgene"
            )

            # get ensids for gene_A
            genesymbols = list(fusion_table['gene_A'])
            ensids = []

            for genesymbol in genesymbols:
                ensid = cv.ensgene2genesymbol(key_file, genesymbol)
                ensids.append(ensid)

            fusion_table['ENSID_A'] = ensids

            # get ensids for gene_B
            genesymbols = list(fusion_table['gene_B'])
            ensids = []

            for genesymbol in genesymbols:
                ensid = cv.ensgene2genesymbol(key_file, genesymbol)
                ensids.append(ensid)

            fusion_table['ENSID_B'] = ensids

            fusion_table.to_csv(
                out_handle,
                sep="\t",
                index=0,
                header=False
            )

        out_handle.close()

        if re.search(r"all|count", args.stage):

            output_file = project_dir + "/" + args.project_name + \
                    "_" + "count_occurrences.txt"
            out_handle = open(output_file, 'w')

            comb_fusions_table = pd.read_csv(
                project_dir + "/" + args.project_name + "_" + "fusion_table.txt",
                sep="\t",
                names=["RANKING",
                   "DETECTION_SOFTWARE",
                   "UNIQ_SAMPLE_ID",
                   "gene_A",
                   "gene_B",
                   "chr_A",
                   "chr_B",
                   "genomicbreakpoint_A",
                   "genomicbreakpoint_B",
                   "splitr_count",
                   "span_count",
                   "score",
                   "FILTER",
                   "UNIQ_SAMPLE_ID_ANNO",
                   "PATIENT_ID",
                   "ENSID_A",
                   "ENSID_B"],
                header=True
            )

            # universe of all variants
            comb_fusions_table = comb_fusions_table.drop_duplicates(
                subset=['gene_A','gene_B',
                        'chr_A', 'chr_B',
                        'genomicbreakpoint_A', 'genomicbreakpoint_B'])

            comb_fusions_table['NUM_OCCURRENCES_IN_COHORT'] = ""
            comb_fusions_table['NAMES_SAMPLES']= ""
            comb_fusions_table['SPLITR_COUNT_COHORT'] = ""
            comb_fusions_table['SPAN_COUNT_COHORT'] = ""

            # and now for each file
            fusion_file = ts.load_tab_delimited(args.input_list)

            for file in fusion_file:
                d={}
                uniq_sample_id = file[0]
                patient_id = file[1]
                file_path = file[2]

                fusion_table = pd.read_csv(
                    file_path,
                    sep="\t",
                    names=["RANKING",
                            "DETECTION_SOFTWARE",
                            "UNIQ_SAMPLE_ID",
                            "gene_A",
                            "gene_B",
                            "chr_A",
                            "chr_B",
                            "genomicbreakpoint_A",
                            "genomicbreakpoint_B",
                            "splitr_count",
                            "span_count",
                            "score",
                            "FILTER"],
                    header=True
                )

                # read in each sample file to nested hash table
                for index, row in fusion_table.iterrows():
                    key_variant = str(row['gene_A']) + ":" + str(row['gene_B']) + \
                        ";" + str(row['chr_A']) + ":" + str(row['chr_B']) + \
                        ";" + str(row['genomicbreakpoint_A']) + ":" + \
                        str(row['genomicbreakpoint_B'])
                    value_variant = {'UNIQ_SAMPLE_ID':row['UNIQ_SAMPLE_ID'],
                                     'splitr_count':row['splitr_count'],
                                     'span_count':row['span_count']}
                    d[key_variant] = value_variant

                # and now fill up the master table
                for index, row in comb_fusions_table.iterrows():
                    key_mut = str(row['gene_A']) + ":" + str(row['gene_B']) + \
                        ";" + str(row['chr_A']) + ":" + str(row['chr_B']) + \
                        ";" + str(row['genomicbreakpoint_A']) + ":" + \
                        str(row['genomicbreakpoint_B'])

                    if key_mut in d:
                        comb_fusions_table.loc[index, 'NAMES_SAMPLES'] = str(d[key_mut]['UNIQ_SAMPLE_ID']) + ";" + str(comb_fusions_table.loc[index, 'NAMES_SAMPLES'])
                        comb_fusions_table.loc[index, 'SPLITR_COUNT_COHORT'] = str(d[key_mut]['splitr_count']) + ";" + str(comb_fusions_table.loc[index, 'SPLITR_COUNT_COHORT'])
                        comb_fusions_table.loc[index, 'SPAN_COUNT_COHORT'] = str(d[key_mut]['span_count']) + ";" + str(comb_fusions_table.loc[index, 'SPAN_COUNT_COHORT'])

            #count the number of occurrences
            for index, row in comb_fusions_table.iterrows():
                comb_fusions_table.loc[index, 'NUM_OCCURRENCES_IN_COHORT'] = len(comb_fusions_table.loc[index,'NAMES_SAMPLES'].split(";")) - 1

            comb_fusions_table = comb_fusions_table[[
                'ENSID_A',
                'ENSID_B',
                'gene_A',
                'gene_B',
                'NUM_OCCURRENCES_IN_COHORT',
                'chr_A',
                'chr_B',
                'genomicbreakpoint_A',
                'genomicbreakpoint_B',
                'NAMES_SAMPLES',
                'SPLITR_COUNT_COHORT',
                'SPAN_COUNT_COHORT'
            ]]

            comb_fusions_table.to_csv(
                out_handle,
                sep="\t",
                index=0,
                header=True
            )

            out_handle.close()
