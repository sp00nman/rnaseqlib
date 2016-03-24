"""
All relevant function to annotate fusions
"""

from rnaseqlib.utils import convert_gene_ids as cv


def map_genes(conversion_table,
              gene_a,
              gene_b):

    a = cv.ensgene2genesymbol(conversion_table, gene_a)
    b = cv.ensgene2genesymbol(conversion_table, gene_b)

    return a, b


def filter_gene_pair(gene_a,
                     gene_b,
                     conversion_table,
                     gene_pairs,
                     tool):
    """
    for each fusion it returns gene A & B
    :param gene_A: upstream gene involved in fusion
    :param gene_B: downstream gene involved in fusion
    :param gene_pairs
    :return:1 present 0 absent
    """

    if tool == "tophatfusion" or tool == "soapfuse":
        (gene_a, gene_b) = map_genes(conversion_table, gene_a, gene_b)

    g = '\t'.join(sorted([gene_a, gene_b]))
    print g

    return "1" if g in gene_pairs else 0