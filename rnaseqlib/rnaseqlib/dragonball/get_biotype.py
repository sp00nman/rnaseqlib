"""
Extract biotype from ensembl database.
"""

import argparse
import logging
import pandas as pd
import gffutils

from rnaseqlib.utils import convert_gene_ids as cv


def get_biotype(gene,
                gtfdb,
                tool,
                conversion_table,
                gene_version):
    """
    Extracts biotype information from database
    :param gene: either up or downstream gene
    :param gtfdb: gtf database
    :param tool: defuse, tophatfusion, soapfuse
    :param conversion_table: converts gene symbol ids to ensembl ids
    :param gene_version: hashtable with version numbers
    :return:
    """

    if tool == "tophatfusion" or tool == "soapfuse":
        gene = cv.ensgene2genesymbol(conversion_table, gene)

    try:
        gene_fromdb = gtfdb[gene + "." + gene_version[gene]]

    except KeyError:
        biotype = "KE"

    else:  # no error occurred
        biotype = "-".join(gene_fromdb.attributes['gene_type'])

    return biotype

