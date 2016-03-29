"""
All relevant function to annotate fusions
"""

from rnaseqlib.utils import convert_gene_ids as cv


def annotate_gene_pair(gene_a,
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
        (gene_a, gene_b) = cv.map_genes(conversion_table, gene_a, gene_b)

    g = '\t'.join(sorted([gene_a, gene_b]))

    return "1" if g in gene_pairs else 0


def annotate_single(gene_a,
                    gene_b,
                    conversion_table,
                    gene_single,
                    tool):
    """

    :param gene_a: upstream gene
    :param gene_b: downstream gene
    :param conversion_table: genesymbol to ensid conversion
    :param gene_single: list of genes
    :param tool: fusion detection tool
    :return:
    """

    if tool == "tophatfusion" or tool == "soapfuse":
        (gene_a, gene_b) = cv.map_genes(conversion_table, gene_a, gene_b)

    if gene_a in gene_single:
        to_return = 'A'

    elif gene_b in gene_single:
        to_return = 'B'

    elif gene_a in gene_single and gene_b in gene_single:
        to_return = 'A,B'
    else:
        to_return = 0

    return to_return