"""
Add column FILTER based on the following criteria
"""


def match_citeria(
        row,
        tool,
        min_num_split_reads,
        min_num_span_reads,
        min_score):
    """

    :param row:
    :param min_num_split_reads:
    :param min_num_span_reads:
    :return:
    """

    filter_annotation = []

    if tool == "defuse" or tool== "tophatfusion":

        if row['score'] < min_score:

            filter_annotation.append("score")

    if int(row['splitr_count']) <= int(min_num_split_reads) or \
        int(row['span_count']) <= int(min_num_span_reads):

        filter_annotation.append("read_evidence")

    if (int(row['readthrough_fusioncatcher']) is 1 or
        int(row['readthrough_ncbi']) is 1 or
        int(row['readthrough_conjoing']) is 1) and \
        (int(row['known_fusion_ticdb']) is 0):

        # the and statement tries to catch fusion like
        # FIP1L1-PDGFRA that are annotated as readthrough fusions,
        # but are mediated through chromosomal deletions

        filter_annotation.append("readthrough")

    if int(row['same_strand_overlapping_fusioncatcher']) is 1 or \
        int(row['same_strand_overlapping_fusioncatcher']) is 1 or \
        int(row['pseudo_genes_fusioncatcher']) is 1 or \
        int(row['partially_overlapping_fusioncatcher']) is 1 or \
        int(row['paralogs_fusioncatcher']) is 1 or \
        int(row['fully_overlapping_fusioncatcher']) is 1:

        filter_annotation.append("false_pos")

    if int(row['healthy_tophatfusion']) is 1 or \
        int(row['healthy_soapfuse']) is 1 or \
        int(row['healthy_fusionmap_less_cons']) is 1 or \
        int(row['healthy_fusionmap']) is 1 or \
        int(row['healthy_fusioncatcher']) is 1 or \
        int(row['healthy_fm_healthy_gene_symbol_less_cons']) is 1 or \
        int(row['healthy_fm_healthy_gene_symbol']) is 1 or \
        int(row['healthy_defuse_c8']) is 1 or \
        int(row['healthy_defuse_c5']) is 1 or \
        int(row['healthy_defuse_a']) is 1:

        # I removed healthy_soapfuse_specific; it's now
        # combined with healthy soapfuse

        filter_annotation.append("healthy")

    if int(row['healthy_validation_pcr']) is 1:

        filter_annotation.append("healthy_pcr")

    if row['3_prime_biotype_ensembl'] != "protein_coding" and \
        row['5_prime_biotype_ensembl'] != "protein_coding":

        filter_annotation.append("no_protein")

    if row['3_prime_biotype_ensembl'] == "pseudogene" and \
        row['3_prime_mapability_encode'] < 0.5:

        filter_annotation.append("biotype_pseudogene")

    if row['5_prime_biotype_ensembl'] == "pseudogene" and \
            row['5_prime_mapability_encode'] < 0.5:

        filter_annotation.append("biotype_pseudogene")

    if row['3_prime_mapability_encode'] < 0.1 or \
        row['5_prime_mapability_encode'] < 0.1:

        filter_annotation.append("mapability")

    if row['5_prime_exclude_encode_duke'] != "NA" or \
        row['5_prime_exclude_encode_dac'] != "NA" or \
        row['3_prime_exclude_encode_duke'] != "NA" or \
        row['3_prime_exclude_encode_dac'] != "NA":

        filter_annotation.append("blacklisted_region")

    if row['distance_ensembl'] < 1000:

        filter_annotation.append("distance_lt1000")

    if tool == "defuse":
        if row['read_through'] == "Y":
            filter_annotation.append("defuse_internal_readthrough")

    if row['5_prime_biotype_ensembl'] is "KE" or \
        row['3_prime_biotype_ensembl'] is "KE":
        filter_annotation.append("deprecated_identifier")

    if row['HLA_imgt|hla_db'] is 'A' or \
        row['HLA_imgt|hla_db'] is 'B' or \
        row['HLA_imgt|hla_db'] is 'A,B':

        filter_annotation.append("HLA_gene")

    if row['HB_ensembl'] is 'A' or \
        row['HB_ensembl'] is 'B' or \
        row['HB_ensembl'] is 'A,B':

        filter_annotation.append("HB_gene")

    if not filter_annotation:
        filter_annotation.append("PASS")

    return ",".join(filter_annotation)
