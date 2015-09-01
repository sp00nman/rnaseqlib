"""
Given a list of genes, extract from table only the genes of interest.
"""

import pandas as pd
from rnaseqlib.utils import tools as ts
from rnaseqlib.utils import convert_gene_ids as cv


def filter_variants(
        df,
        ExonicFunc_ensGene="nonsynonymous_SNV",
        QUAL=0,
        DP=0,
        QD=0
):
    """
    Simple function to filter for low confidence variants.
    :param input_file: Input file
    :param ExonicFunc_ensGene: Ensembl gene ids
    :param QUAL: Variant quality
    :param DP: Variant read depth; number of reads covering that position
    :param QD: Variant quality by depth
    :return: pandas dataframe object
    """

    df = df[df['ExonicFunc_ensGene'] == ExonicFunc_ensGene]
    df = df[df['QUAL'] > QUAL]
    df = df[df['DP'] > DP]
    df = df[df['QD'] > QD]

    return df


def extract_genes(project_dir,
                  output_file,
                  target_list,
                  vcf2table,
                  conversion_table):
    """

    :param project_dir: Project directory
    :param output_file: Name of output file
    :param target_list: List of genes of interest [ENSID\nENSID..]
    :param vcf2table: Dataframe with selected columns.
    :return:None
    """

    vcfs = ts.load_tab_delimited(vcf2table)
    targets = [line.rstrip('\n') for line in open(target_list, 'r')]
    out_handle = open(output_file, 'a')

    for vcf_file in vcfs:

        uniq_sample_id = vcf_file[0]
        patient_id = vcf_file[1]

        print uniq_sample_id
        print patient_id

        print project_dir + "_" + uniq_sample_id + ".vcf2table"

        variant_table = pd.read_csv(
            project_dir + "_" + uniq_sample_id + ".vcf2table",
            sep="\t",
            names=["CHROM", "POS", "REF", "ALT", "QUAL",
                   "DP", "QD", "Gene_ensGene", "ExonicFunc_ensGene",
                   "GT", "AD"],
            header=True
        )

        #filtered_table = filter_variants(
        #    df=variant_table,
        #    ExonicFunc_ensGene="non_synonymous_SNV",
        #    QUAL=0,
        #    DP=0,
        #    QD=0
        #)

        variant_table = variant_table[variant_table['ExonicFunc_ensGene'] == "nonsynonymous_SNV"]
        variant_table = variant_table[variant_table['QUAL'] > 50]
        variant_table = variant_table[variant_table['DP'] > 5]
        variant_table = variant_table[variant_table['QD'] > 5]


        selected_genes = variant_table[variant_table.Gene_ensGene.isin(
            targets)]

        # add patient information
        selected_genes['UNIQ_SAMPLE_ID'] = uniq_sample_id
        selected_genes['PATIENT_ID'] = patient_id

        # append to file
        selected_col = selected_genes[['UNIQ_SAMPLE_ID',
                                       'PATIENT_ID',
                                       'Gene_ensGene',
                                       'AD']]
        # convert gene ids
        key_file = cv.read_ensgene_genesymb(conversion_table)
        ensids = list(selected_col['Gene_ensGene'])
        genesymbols = []
        for ensid in ensids:
            genesymbol = cv.ensgene2genesymbol(key_file, ensid)
            genesymbols.append(genesymbol)

        # add column
        selected_col['GENESYMBOL'] = genesymbols

        selected_col.to_csv(
            out_handle,
            sep="\t",
            index=0,
            header=False
        )

    out_handle.close()
