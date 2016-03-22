"""
All relevant function to annotate fusions
"""


def get_map_genes(refseqid_ensemblid, ensgene_genesymb,
                  line, options, col_num_gene1, col_num_gene2):

    if options.tool == 'defuse':
        a = line[col_num_gene1[options.tool]]
        b = line[col_num_gene2[options.tool]]

    elif options.tool == 'fusionmap':
        print line
        a = loop_fusionmap_refseqids(refseqid_ensemblid, line,
                                     options, col_num_gene1)
        b = loop_fusionmap_refseqids(refseqid_ensemblid, line,
                                     options, col_num_gene2)

    else:  # tophat-fusion
        a = genesymb_to_ensgene(ensgene_genesymb,
                                line[col_num_gene1[options.tool]])
        b = genesymb_to_ensgene(ensgene_genesymb,
                                line[col_num_gene2[options.tool]])

    return a, b



def filter_gene_pair(gene_A,
                     gene_B,
                     gene_pairs):
    """
    for each fusion it returns gene A & B
    :param gene_A:
    :param gene_B:
    :param gene_pairs
    :return:1 present 0 absent
    """

    (a, b) = get_map_genes(refseqid_ensemblid,
                           ensemblid,
                           gene_A,
                           gene_B)

    g = '\t'.join(sorted([a, b]))

    return "1" if g in gene_pairs else 0