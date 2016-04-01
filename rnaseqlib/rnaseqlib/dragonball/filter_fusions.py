"""
Add column FILTER based on the following criteria
"""


def match_citeria(row):
    """

    :param fusions: dataframe with PASS or
    :return: str; with the filter annotation
    """

    filter_annotation = ''

    if row['score'] < 0.8:

        filter_annotation = "score"

    if row['readthrough_fusioncatcher'] == 1 or \
        row['readthrough_ncbi'] == 1 or \
        row['readthrough_conjoing'] == 1:

        filter_annotation = filter_annotation + "," + "readthrough"

    if row['same_strand_overlapping_fusioncatcher'] == 1 or \
        row['same_strand_overlapping_fusioncatcher'] == 1 or \
        row['pseudo_genes_fusioncatcher'] == 1 or \
        row['partially_overlapping_fusioncatcher'] == 1 or \
        row['paralogs_fusioncatcher'] == 1 or \
        row['fully_overlapping_fusioncatcher'] == 1:

        filter_annotation = filter_annotation + "," + "false_pos"

    if row['healthy_tophatfusion'] == 1 | \
        row['healthy_soapfuse_specific'] == 1 | \
        row['healthy_soapfuse'] == 1 | \
        row['healthy_fusionmap_less_cons'] == 1 | \
        row['healthy_fusionmap'] == 1 | \
        row['healthy_fusioncatcher'] == 1 | \
        row['healthy_fm_healthy_gene_symbol_less_cons'] == 1 | \
        row['healthy_fm_healthy_gene_symbol'] == 1 | \
        row['healthy_defuse_c8'] == 1 | \
        row['healthy_defuse_c5'] == 1 | \
        row['healthy_defuse_a'] == 1:

        filter_annotation = filter_annotation + "," + "healthy"

    if row['3_prime_biotype_ensembl'] != "protein_coding" and \
        row['5_prime_biotype_ensembl'] != "protein_coding":

        filter_annotation = filter_annotation + "," + "no_protein"

    if row['3_prime_mapability_encode'] < 0.1 or \
        row['5_prime_mapability_encode'] < 0.1:

        filter_annotation = filter_annotation + "," + "mapability"

    if row['5_prime_exclude_encode_duke'] != "NA" or \
        row['5_prime_exclude_encode_dac'] != "NA" or \
        row['3_prime_exclude_encode_duke'] != "NA" or \
        row['3_prime_exclude_encode_dac'] != "NA":

        filter_annotation = filter_annotation + "," + "blacklisted_region"

    if row['distance_ensembl'] < 1000:

        filter_annotation = filter_annotation + "," + "distance_lt1000"

    if filter_annotation == '':
        filter_annotation = "PASS"

    return filter_annotation
