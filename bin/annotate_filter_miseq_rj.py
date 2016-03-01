#!/usr/bin/env python

import argparse
import re
import os
import logging
import pandas as pd

from rnaseqlib.utils import tools as ts
from rnaseqlib.varcall import annovar_runnables as annovar_rb
from rnaseqlib.varcall import gatk_runnables as gatk
from rnaseqlib.varcall import filter_vcf as fv
from rnaseqlib.utils import convert_gene_ids as cv
from rnaseqlib.varcall import extract_genes_of_interest as goi


VCF_TABLE = \
            ["CHROM",
             "POS",
             "ID",
             "REF",
             "ALT",
             "SNPEFF_GENE_NAME",
             "Gene.ensGene",
             "SNPEFF_GENE_BIOTYPE",
             "Func.ensGen",
             "ExonicFunc.ensGene",
             "AAChange.ensGene",
             "QUAL",
             "QD",
             "VQSLOD",
             "cytoBand",
             "genomicSuperDups",
             "snp129",
             "snp142Common",
             "1000g2015feb_all",
             "1000g2015feb_eur",
             "1000g2015feb_afr",
             "1000g2015feb_amr",
             "1000g2015feb_eas",
             "1000g2015feb_sas",
             "esp5400_all",
             "esp6500siv2_all",
             "cosmic70",
             "clinvar_20150330",
             "SIFT_score",
             "SIFT_pred",
             "Polyphen2_HDIV_score",
             "Polyphen2_HDIV_pred",
             "Polyphen2_HVAR_score",
             "Polyphen2_HVAR_pred",
             "LRT_score",
             "LRT_pred",
             "MutationTaster_score",
             "MutationTaster_pred",
             "FATHMM_score",
             "FATHMM_pred",
             "RadialSVM_score",
             "RadialSVM_pred",
             "LR_score",
             "LR_pred",
             "VEST3_score",
             "CADD_raw",
             "CADD_phred",
             "GERP++_RS",
             "phyloP46way_placental",
             "phyloP100way_vertebrate",
             "SiPhy_29way_logOdds",
             "caddgt10",
             "CADDindel",
             "CADDindel_Phred",
             "GT",
             "AD",
             "DP",
             "GQ",
             "PL"]

SELECT_COLUMNS = \
            ["UNIQ_SAMPLE_ID",
             "CHROM",
             "POS",
             "ID",
             "REF",
             "ALT",
             "GT",
             "AD",
             "DP",
             "GQ",
             "PL",
             "GENESYMBOL",
             "ENSEMBL_GENEID",
             "ENSEMBL_TRANSCRIPTID",
             "Func.ensGen",
             "ExonicFunc.ensGene",
             "QUAL",
             "QD",
             "CANONICAL_TRANSCRIPT_EXON_NUM",
             "CANONICAL_TRANSCRIPT_NUCLEOTIDE_CHANGE",
             "CANONICAL_TRANSCRIPT_AS_CHANGE",
             "cytoBand",
             "genomicSuperDups",
             "snp129",
             "snp142Common",
             "1000g2015feb_all",
             "1000g2015feb_eur",
             "1000g2015feb_afr",
             "1000g2015feb_amr",
             "1000g2015feb_eas",
             "1000g2015feb_sas",
             "esp5400_all",
             "esp6500siv2_all",
             "cosmic70",
             "clinvar_20150330",
             "SIFT_score",
             "SIFT_pred",
             "Polyphen2_HDIV_score",
             "Polyphen2_HDIV_pred",
             "Polyphen2_HVAR_score",
             "Polyphen2_HVAR_pred",
             "LRT_score",
             "LRT_pred",
             "MutationTaster_score",
             "MutationTaster_pred",
             "FATHMM_score",
             "FATHMM_pred",
             "RadialSVM_score",
             "RadialSVM_pred",
             "LR_score",
             "LR_pred",
             "VEST3_score",
             "CADD_raw",
             "CADD_phred",
             "GERP++_RS",
             "phyloP46way_placental",
             "phyloP100way_vertebrate",
             "SiPhy_29way_logOdds",
             "caddgt10",
             "CADDindel",
             "CADDindel_Phred",
            ]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Annotate vcf files for RJ.')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'analysis stage.[process_all,annovar,selection,'
                             'clonalityXchr, vcf2table_rj,'
                             'vcf2table_dna, reformat_rj]')
    parser.add_argument('--project_name', required=False, type=str,
                        help="name of the project")
    parser.add_argument('--patient_vcf', required=False, type=str,
                        help='List of VCF files/patient to process.')

    parser.add_argument('--sample_dir', required=False, type=str,
                        help="Path to sample directory")
    parser.add_argument('--output_dir', required=False, type=str,
                        help="Path to output directory.")
    # resources
    parser.add_argument('--ref_genome', required=False, type=str,
                        help='Reference genome.')
    parser.add_argument('--id_conversion', required=False, type=str,
                        help='Conversion table. Ensembl ids to genesymbols')
    parser.add_argument('--canonical_transcripts', required=False, type=str,
                        help='Select canonical transcript.')
    parser.add_argument('--mae_genes', required=False, type=str,
                        help='Table of monoallelic expressed genes based on'
                             'a study by Savova et al., Nature Genetics, 2015')

    # annovar specific options
    parser.add_argument('--annovar', required=False, type=str,
                        help="Annotate variant with annovar.")

    # hardware specific options
    parser.add_argument('--num_cpus', dest='num_cpus', required=False,
                        help='Number of cpus.')
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
    logging.info("Annotate and select vcf files.")

    # start workflow
    if re.search(r"process_all|annovar", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            cmd = annovar_rb.run_annovar(
                vcf_file=path2vcf,
                protocol="ensGene,"
                        + "cytoBand,"
                        + "genomicSuperDups,"
                        + "snp142Mult,"
                        + "snp129,"
                        + "snp142Common,"
                        + "1000g2015feb_all,"
                        + "1000g2015feb_afr,"
                        + "1000g2015feb_amr,"
                        + "1000g2015feb_eas,"
                        + "1000g2015feb_eur,"
                        + "1000g2015feb_sas,"
                        + "esp5400_all,"
                        + "esp6500siv2_all,"
                        + "cosmic70,"
                        + "clinvar_20150330,"
                        + "ljb26_all,"
                        + "caddgt10,"
                        + "caddindel",
                output_file=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id,
                operation="g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f",
                buildversion="hg19",
                annovar_dir=args.annovar
            )

            status = ts.run_cmd(
                message=stdout_msg['annotation'],
                command=cmd,
                debug=args.debug
            )

    if re.search(r"process_all|selection", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            fv.filter_vcf(
                input_file=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id + "."
                        + "hg19_multianno.vcf",
                file_prefix=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id,
                file_suffix=file_ext['inhouse'],
                mode="SOMATIC"
            )

    if re.search(r"process_all|clonalityXchr", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            fv.filter_vcf(
                input_file=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id + "."
                        + "hg19_multianno.vcf",
                file_prefix=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id,
                file_suffix=file_ext['inhouse'],
                mode="CLONXCHR"
            )

    if re.search(r"process_all|mae_autosomes", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            print uniq_sample_id

            fv.filter_vcf(
                input_file=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id + "."
                        + "hg19_multianno.vcf",
                file_prefix=project_dir + "/"
                        + args.project_name + "."
                        + uniq_sample_id,
                file_suffix=file_ext['inhouse'],
                mode="MAE",
                ensids=ts.load_mae_genes(args.mae_genes)
            )

    if re.search(r"process_all|vcf2table_rj", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            cmd = gatk.var2table_proudpv(
                input_file=project_dir + "/"
                           + args.project_name + "."
                           + uniq_sample_id + "."
                           + file_ext['proud_pv'],
                output_file=project_dir + "/" + args.project_name + "."
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

    if re.search(r"process_all|vcf2table_dna", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]
            path2vcf = vcf_file[2]

            # TODO: automatic change from proud_pv to clonewars...
            #input_file=project_dir + "/"
            #    + args.project_name + "."
            #    + uniq_sample_id + "."
            #    + file_ext['clonewars']

            # normally
            # input_file=project_dir + "/"
            #+ args.project_name + "."
            # + uniq_sample_id + "."
            #+ file_ext['proud_pv'],

            cmd = gatk.var2table_dna(
                input_file=project_dir + "/"
                           + args.project_name + "."
                           + uniq_sample_id + "."
                           + file_ext['mae'],
                output_file=project_dir + "/" + args.project_name + "."
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

    if re.search(r"process_all|reformatXchr", args.stage):

        vcfs = ts.load_tab_delimited(args.patient_vcf)
        canonical_transcripts = ts.load_dictionary(
            args.canonical_transcripts,
            sep="\t")

        for vcf_file in vcfs:

            uniq_sample_id = vcf_file[0]
            patient_id = vcf_file[1]

            # set header to index 0 or 1
            variant_table = pd.read_csv(
                project_dir + "/" + args.project_name + "."
                + uniq_sample_id + ".vcf2table.txt",
                sep="\t",
                names=VCF_TABLE,
                header=True
            )

            # add patient information
            variant_table['UNIQ_SAMPLE_ID'] = uniq_sample_id

            # convert gene ids
            key_file = cv.read_ensgene_genesymb(args.id_conversion)
            ensids = list(variant_table['Gene.ensGene'])
            genesymbols = []
            for ensid in ensids:
                genesymbol = cv.ensgene2genesymbol(key_file, ensid)
                genesymbols.append(genesymbol)

            # add column
            variant_table['GENESYMBOL'] = genesymbols

            # only select canonical transcript
            variant_table['CANONICAL_TRANSCRIPT_AAChange'] = \
                variant_table.apply(lambda row: goi.select_canonical_transcript(
                    row['AAChange.ensGene'],
                    canonical_transcripts=canonical_transcripts,
                    sep=","), axis=1
                )

            # split AAChange_ensGene by ':'
            variant_table_split = variant_table[
                'CANONICAL_TRANSCRIPT_AAChange'].apply(lambda row: pd.Series(row.split(':')))
            # rename columns
            variant_table_split.columns = [
                'ENSEMBL_GENEID',
                'ENSEMBL_TRANSCRIPTID',
                'CANONICAL_TRANSCRIPT_EXON_NUM',
                'CANONICAL_TRANSCRIPT_NUCLEOTIDE_CHANGE',
                'CANONICAL_TRANSCRIPT_AS_CHANGE']

            # concatenate the dataframes
            variant_table_concat = pd.concat(
                [variant_table,variant_table_split], axis=1)

            # select specific columns
            variant_table_concat = variant_table_concat[SELECT_COLUMNS]

            output_file = project_dir + "/" + args.project_name + "." \
                          + uniq_sample_id + ".variant_table.txt"

            out_handle = open(output_file, 'w')
            variant_table_concat.to_csv(
                out_handle,
                sep="\t",
                index=0,
                header=True
            )
            out_handle.close()

    if re.search(r"process_all|reformat_rj", args.stage):

        canonical_transcripts = ts.load_dictionary(
            args.canonical_transcripts,
            sep="\t")

        vcfs = ts.load_tab_delimited(args.patient_vcf)

        for vcf_file in vcfs:

            print vcf_file
            uniq_sample_id = vcf_file[0]

            # set header to index 0 or 1
            variant_table = pd.read_csv(
                project_dir + "/" + args.project_name + "."
                + uniq_sample_id + ".vcf2table.txt",
                sep="\t"
            )

            # add patient information
            variant_table['UNIQ_SAMPLE_ID'] = uniq_sample_id

            # convert ensembl gene ids to genesymbols
            key_file = cv.read_ensgene_genesymb(args.id_conversion)
            ensids = list(variant_table['Gene.ensGene'])
            genesymbols = []

            for ensid in ensids:
                genesymbol = cv.ensgene2genesymbol(key_file, ensid)
                genesymbols.append(genesymbol)

            variant_table['GENESYMBOL'] = genesymbols

            # select canonical transcript
            if not variant_table.empty:

                canonical_transcripts = \
                    variant_table.apply(
                        lambda row: pd.Series(goi.select_canonical_transcript(
                        row['AAChange.ensGene'],
                        canonical_transcripts=canonical_transcripts,
                        sep=",")), axis=1
                    )
                canonical_transcripts = canonical_transcripts.rename(
                    columns={0:'CANONICAL_TRANSCRIPT_AAChange',
                             1:'IS_CANONICAL_TRANSCRIPT'})
                variant_table_merge = pd.concat(
                    [variant_table, canonical_transcripts], axis=1)

                # split AAChange_ensGene by ':'
                variant_table_split = variant_table_merge[
                    'CANONICAL_TRANSCRIPT_AAChange'].apply(
                    lambda row: pd.Series(row.split(':')))

                # rename columns
                variant_table_split.columns = [
                    'ENSEMBL_GENEID',
                    'ENSEMBL_TRANSCRIPTID',
                    'CANONICAL_TRANSCRIPT_EXON_NUM',
                    'CANONICAL_TRANSCRIPT_NUCLEOTIDE_CHANGE',
                    'CANONICAL_TRANSCRIPT_AS_CHANGE']

                # concatenate the dataframes
                variant_table_concat = pd.concat(
                    [variant_table_merge, variant_table_split], axis=1)

            output_file = project_dir + "/" \
                          + args.project_name + "." \
                          + uniq_sample_id \
                          + ".variant_table.txt"

            out_handle = open(output_file, 'w')

            # TODO: variable can not be defined if variant_table is empty
            variant_table_concat.to_csv(
                out_handle,
                sep="\t",
                index=0,
                header=True
            )
            out_handle.close()