"""
Calculates the distance between genes on the fly.
"""

import argparse
import logging
import pandas as pd
import gffutils

from rnaseqlib.utils import convert_gene_ids as cv


def min_dis_gene(gene_a,
                 gene_b,
                 conversion_table,
                 tool,
                 gtfdb,
                 gene_version):
    """

    :param gene_a:
    :param gene_b:
    :param ens2gs_conversion_table:
    :param annotation_dta:
    :param tool:
    :return:
    """

    if tool == "tophatfusion" or tool == "soapfuse":
        (gene_a, gene_b) = cv.map_genes(conversion_table, gene_a, gene_b)

    try:
        up_gene = gtfdb[gene_a + "." + gene_version[gene_a]]
        down_gene = gtfdb[gene_b + "." + gene_version[gene_b]]

    except KeyError:
        min_dis = "KE"

    else:  # no error occurred

        if up_gene.chrom == down_gene.chrom and up_gene.strand == down_gene.strand:
            min_dis = min(
                abs(up_gene.start - down_gene.start),
                abs(up_gene.start - down_gene.end),
                abs(up_gene.end - down_gene.start),
                abs(up_gene.end - down_gene.end))
        else:
            min_dis = "NA"

    return min_dis


