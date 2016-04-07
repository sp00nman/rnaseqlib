"""
Collection of functions that execute external software or UNIX commands.
"""


def bam_stat(input_file,
             output_file,
             map_qual=30):
    """
    Calculates reads mapping statistics for a given BAM file
    http://rseqc.sourceforge.net/#bam_stat.py
    :param input_file: BAM file
    :param map_qual: minimum mapping quality for an alignment to be called
    "uniquely mapped". default=30
    :return:str object; Command to be executed
    """

    cmd_bam_stat = "bam_stat.py"
    cmd_bam_stat += " -i " + input_file
    cmd_bam_stat += " -q " + str(map_qual)
    cmd_bam_stat += " >>" + output_file
    cmd_bam_stat += " 2>&1"  # redirects stdout and stderror to output file

    return cmd_bam_stat


def genebody_coverage(input_file,
                      reference,
                      output_prefix):
    """
    Read coverage over gene body. This module is used to check if reads
    coverage is uniform and if there is any 5p or 3p bias.
    Default values as described http://rseqc.sourceforge.net/#genebody_coverage.py
    :param input_file: BAM file
    :param reference: Reference gene model in BED format
    :param output_prefix: Prefix of output file
    :return: str object; Command to be executed
    """

    cmd_genebody_coverage = "geneBody_coverage.py"
    cmd_genebody_coverage += " -i " + input_file
    cmd_genebody_coverage += " -r " + reference
    cmd_genebody_coverage += " -o " + output_prefix

    return cmd_genebody_coverage


def inner_distance(input_file,
                   reference,
                   output_prefix):
    """
    Calculate the inner distance (or insert size) between two paired RNA reads.
    Default values as described: http://rseqc.sourceforge.net/#inner_distance.py
    :param input_file: BAM file
    :param reference: Reference gene model in BED format.
    :param output_prefix: Prefix of output file
    :return: str object; Command to be executed
    """

    cmd_inner_distance = "inner_distance.py"
    cmd_inner_distance += " -i " + input_file
    cmd_inner_distance += " -r " + reference
    cmd_inner_distance += " -o " + output_prefix

    return cmd_inner_distance


def junction_annotation(input_file,
                        reference,
                        output_prefix):
    """
    Compares detected splice junctions to reference gene model
    Default values as described: http://rseqc.sourceforge.net/#junction_annotation.py
    :param input_file: BAM file
    :param reference: Reference gene model in BED format.
    :param output_prefix: Prefix of output file
    :return: str object; Command to be executed
    """

    cmd_junction_annotation = "junction_annotation.py"
    cmd_junction_annotation += " -i " + input_file
    cmd_junction_annotation += " -r " + reference
    cmd_junction_annotation += " -o " + output_prefix

    return cmd_junction_annotation


def junction_saturation(input_file,
                        reference,
                        output_prefix):
    """
    Checks for saturation by resampling 5%, 10%, 15%, ..., 95% of total alignments
    Default values as described: http://rseqc.sourceforge.net/#junction_saturation.py
    :param input_file: BAM file
    :param reference: Reference gene model in BED format.
    :param output_prefix: Prefix of output file
    :return: str object; Command to be executed
    """

    cmd_junction_saturation = "junction_saturation.py"
    cmd_junction_saturation += " -i " + input_file
    cmd_junction_saturation += " -r " + reference
    cmd_junction_saturation += " -o " + output_prefix

    return cmd_junction_saturation


def read_duplication(input_file,
                     output_prefix):
    """
    Calculates duplication rates
    Default values as described: http://rseqc.sourceforge.net/#read_duplication.py
    :param input_file: BAM file
    :param output_prefix: Prefix of output file
    :return: str object; Command to be executed
    """

    cmd_read_duplication = "read_duplication.py"
    cmd_read_duplication += " -i " + input_file
    cmd_read_duplication += " -o " + output_prefix

    return cmd_read_duplication


def read_gc(input_file,
            output_prefix):
    """
    Calculates GC content
    Default values as described: http://rseqc.sourceforge.net/#read_GC.py
    :param input_file: BAM file
    :param output_prefix: Prefix of output file
    :return: str object; Command to be executed
    """

    cmd_read_gc = "read_GC.py"
    cmd_read_gc += " -i " + input_file
    cmd_read_gc += " -o " + output_prefix

    return cmd_read_gc