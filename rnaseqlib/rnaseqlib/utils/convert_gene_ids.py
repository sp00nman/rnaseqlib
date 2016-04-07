"""
Helper function to map ensembl gene ids to genesymbols and vice versa.
ensgene - ensembl gene id (eg. ENSG....)
genesymb - gene symbols (eg. JAK2)
"""

from rnaseqlib.utils import tools as ts


def read_ensgene_genesymb(filename,
                          mode="ensgene2genesymbol"):
    """
    *file structure*
    ensgene_id      gene_symbol
    ENSG00000000003 TSPAN6

    *for mode=="ensgene2genesymbol"*
    key|value key=ensgeneid | value=gene symbol

    *for mode=="genesymbol2ensgene"

    :param filename: name of filename with matching ids
    :param mode: [ensgene2genesymbol, genesymbol2ensgene]
    :return:dictionary with
    """
    ensgene_genesymb = ts.load_tab_delimited(filename)
    tab_ensgene_genesymb = {}

    if mode == "ensgene2genesymbol":
        for line in ensgene_genesymb:
            ensgeneid = line[0]
            genesymb = line[1]
            tab_ensgene_genesymb[ensgeneid] = genesymb

    if mode == "genesymbol2ensgene":
        for line in ensgene_genesymb:
            ensgeneid = line[0]
            genesymb = line[1]
            tab_ensgene_genesymb[genesymb] = ensgeneid

    return tab_ensgene_genesymb


def ensgene2genesymbol(tab_ensgene_genesymb,
                       ensgene):
    """
    Matches either genesymbols or ensids to its corresponding partner
    careful: LOC541471 changed ?
    """
    if ensgene in tab_ensgene_genesymb:
        ensgeneid = tab_ensgene_genesymb[ensgene]

    else:
        # keep gene_symbol if nothing was found
        ensgeneid = ensgene
    return ensgeneid


def map_genes(conversion_table,
              gene_a,
              gene_b):
    """
    Questionable weather this function is necessary, seams like
    an useless wrapper
    :param conversion_table: genesymbol to ensid conversion
    :param gene_a: upstream gene
    :param gene_b: downstream gene
    :return:
    """

    a = ensgene2genesymbol(
        conversion_table, gene_a)
    b = ensgene2genesymbol(
        conversion_table, gene_b)

    return a, b