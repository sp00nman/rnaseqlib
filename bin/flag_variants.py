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

    parser = argparse.ArgumentParser(description='Flag variants.')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'analysis stage.[all]')
    parser.add_argument('--output_dir', required=False, type=str,
                        help='Project directory output.')
    parser.add_argument('--project_name', required=False, type=str,
                        help='Name of the project.')
    parser.add_argument('--vcf2flag', required=False, type=str,
                        help='List of VCF files used for flagging')

    # resources
    parser.add_argument('--ref_genome', required=False, type=str,
                        help='Reference genome.')
    parser.add_argument('--id_conversion', required=False, type=str,
                        help='Conversion table. Ensembl ids to genesymbols')

    # java related settings
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
    logging.info("Flag variants")

    # start workflow
    if re.search(r"all|vcf2table", args.stage):

        vcfs = ts.load_tab_delimited(args.vcf2flag)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            cmd = gatk.variants2table(
                input_file=path2vcf,
                output_file=project_dir + "/" + args.project_name + "_"
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

    if re.search(r"all|patient_name2col", args.stage):

        output_file = project_dir + "/" + args.project_name + "_" \
                      + "variant_table.txt"
        vcfs = ts.load_tab_delimited(args.vcf2flag)
        out_handle = open(output_file, 'a')

        for vcf_file in vcfs:
            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]

            variant_table = pd.read_csv(
            project_dir + "/" + args.project_name + "_" \
            + uniq_sample_id + ".vcf2table",
            sep="\t",
            names=["CHROM", "POS", "REF", "ALT", "QUAL",
                   "DP", "QD", "Gene_ensGene", "ExonicFunc_ensGene",
                   "GT", "AD"],
            header=True
            )

            # add patient information
            variant_table['UNIQ_SAMPLE_ID'] = uniq_sample_id
            variant_table['PATIENT_ID'] = patient_id

            # convert gene ids
            key_file = cv.read_ensgene_genesymb(args.id_conversion)
            ensids = list(variant_table['Gene_ensGene'])
            genesymbols = []
            for ensid in ensids:
                genesymbol = cv.ensgene2genesymbol(key_file, ensid)
                genesymbols.append(genesymbol)

            # add column
            variant_table['GENESYMBOL'] = genesymbols

            variant_table.to_csv(
                out_handle,
                sep="\t",
                index=0,
                header=False
            )

        out_handle.close()

    if re.search(r"all|count", args.stage):

        input_file = project_dir + "/" + args.project_name + \
                    "_" + "count_occurences.txt"
        out_handle = open(input_file, 'w')

        mutation_table = pd.read_csv(
            project_dir + "/" + args.project_name + "_" + "variant_table.txt",
            sep="\t",
            names=["CHROM", "POS", "REF", "ALT", "QUAL",
                   "DP", "QD", "Gene_ensGene", "ExonicFunc_ensGene",
                   "GT", "AD", "UNIQ_SAMPLE_ID", "PATIENT_ID", "GENESYMBOL"],
            header=True
        )

        # universe of all variants
        mutation_table = mutation_table.drop_duplicates(
            subset=['CHROM','POS','REF', 'ALT'])

        mutation_table['DP_COHORT'] = ""
        mutation_table['QUAL_COHORT'] = ""
        mutation_table['QD_COHORT'] = ""
        mutation_table['Gene_ensGene_COHORT'] = ""
        mutation_table['ExonicFunc_ensGene_COHORT'] = ""
        mutation_table['GT_COHORT'] = ""
        mutation_table['AD_COHORT'] = ""
        mutation_table['UNIQ_SAMPLE_ID_COHORT'] = ""

        # and now for each vcf
        vcfs = ts.load_tab_delimited(args.vcf2flag)

        for vcf_file in vcfs:
            d = {} # for each vcf file
            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]

            print uniq_sample_id
            variant_table = pd.read_csv(
                project_dir + "/" + args.project_name + "_" \
                + uniq_sample_id + ".vcf2table",
                sep="\t",
                names=["CHROM", "POS", "REF", "ALT", "QUAL",
                       "DP", "QD", "Gene_ensGene", "ExonicFunc_ensGene",
                       "GT", "AD"],
                header=True
            )

            # read in vcf to nested hash table
            for index, row in variant_table.iterrows():
                key_variant = str(row['CHROM']) + ":" + str(row['POS']) + ";" \
                      + row['REF'] + ";" + row['ALT']
                value_variant = {'QUAL':row['QUAL'],
                         'DP':row['DP'],
                         'QD':row['QD'],
                         'Gene_ensGene':row['Gene_ensGene'],
                         'ExonicFunc_ensGene':row['ExonicFunc_ensGene'],
                         'GT':row['GT'],
                         'AD':row['AD']}
                d[key_variant] = value_variant

            for index, row in mutation_table.iterrows():
                key_mut = str(row['CHROM']) + ":" + str(row['POS']) + ";" \
                          + row['REF'] + ";" + row['ALT']

                if key_mut in d:
                    #print d[key_mut]
                    #print d[key_mut]['DP']
                    mutation_table.loc[index,'QUAL_COHORT'] = str(d[key_mut]['QUAL']) + ";" + str(mutation_table.loc[index, 'QUAL_COHORT'])
                    mutation_table.loc[index,'DP_COHORT'] = str(d[key_mut]['DP']) + ";" + str(mutation_table.loc[index,'DP_COHORT'])
                    mutation_table.loc[index,'QD_COHORT'] = str(d[key_mut]['QD']) + ";" + str(mutation_table.loc[index,'QD_COHORT'])
                    mutation_table.loc[index,'ExonicFunc_ensGene_COHORT'] = str(d[key_mut]['ExonicFunc_ensGene']) + ";" + str(mutation_table.loc[index,'ExonicFunc_ensGene_COHORT'])
                    mutation_table.loc[index,'GT_COHORT'] = str(d[key_mut]['GT']) + ";" + str(mutation_table.loc[index,'GT_COHORT'])
                    mutation_table.loc[index,'AD_COHORT'] = str(d[key_mut]['AD']) + ";" + str(mutation_table.loc[index,'AD_COHORT'])
                    mutation_table.loc[index,'UNIQ_SAMPLE_ID_COHORT'] = uniq_sample_id + ";" + str(mutation_table.loc[index,'UNIQ_SAMPLE_ID_COHORT'])

        # select columns
        mutation_table = mutation_table[[
            'CHROM',
            'POS',
            'REF',
            'ALT',
            'Gene_ensGene',
            'GENESYMBOL',
            'UNIQ_SAMPLE_ID_COHORT',
            'QUAL_COHORT',
            'DP_COHORT',
            'QD_COHORT',
            'ExonicFunc_ensGene_COHORT',
            'GT_COHORT',
            'AD_COHORT'
        ]]

        mutation_table.to_csv(
            out_handle,
            sep="\t",
            index=0,
            header=False
        )

        out_handle.close()






