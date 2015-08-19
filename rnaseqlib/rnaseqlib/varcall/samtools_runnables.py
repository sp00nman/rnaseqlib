"""
Collection of functions that execute external software or UNIX commands.
"""


def extract(input_file,
            region,
            output_file):
    """
    Extract a region from a BAM file.
    """

    cmd_extract = "samtools view -b -h %s %s >%s " % (input_file,
                                                      region,
                                                      output_file)
    return cmd_extract


def index_bam(inbamfile):
    """
    Index BAM file.
    :param inbamfile: name of BAM formatted file
    :return: Command to be executed; type str
    """

    cmd_index_bam = "samtools index %s" % (inbamfile)

    return cmd_index_bam


def varcall_samtools(input_file,
                     output_file,
                     ref_genome):
    """
    Variant calling with samtools.
    """

    cmd_samtools = "samtools mpileup " \
                   "-C 50 " \
                   "-uf %s " \
                   "%s " \
                   "| bcftools view -vcg - > %s" % (ref_genome,
                                                    input_file,
                                                    output_file)
    return cmd_samtools
