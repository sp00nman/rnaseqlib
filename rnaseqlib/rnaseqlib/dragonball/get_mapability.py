"""
For each breakpoint (or position X) calculate the mean mapability of
the region [x-50 to x+50]
"""

import pybedtools


def get_mean_values(
        bigwigfile,
        chr,
        breakpoint,
        tool):
    """

    :param bigwigfile: file in BIGWIG format
    :param chr: chromosome of gene X
    :param breakpoint: genomic breakpoint of gene X
    :return: str? mapability score X for region [x-50 to x+50]
    """

    if tool == "defuse":
        chromosome = "chr" + str(chr)

    else:
        chromosome = str(chr)

    mean_score = bigwigfile.stats(
        str(chromosome),
        breakpoint - 50,
        breakpoint + 50)

    return mean_score[0]


def intersect_genomic_breakpoint(
        bedfile,
        chr,
        breakpoint,
        tool):
    """

    :param bedfile: file in BED format
    :param chr: chromosome of gene X
    :param breakpoint: genomic breakpoint of gene X
    :return: str with intersected region
    """
    if tool == "defuse":
        query = "chr" + str(chr) + " " + str(breakpoint) + " " + str(breakpoint + 1)
    else:
        query = str(chr) + " " + str(breakpoint) + " " + str(breakpoint + 1)

    interval = pybedtools.BedTool(query, from_string=True)
    result = bedfile.intersect(interval)

    if result:
        return_value = str(result[0].chrom) + ":" \
                       + str(result[0].start) + "-" \
                       + str(result[0].end) + ";" \
                       + str(result[0].name)
    else:
        return_value = "NA"

    return return_value