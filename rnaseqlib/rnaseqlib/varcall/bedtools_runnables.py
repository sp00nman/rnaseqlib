"""
Collection of functions that execute external software or UNIX commands.
"""


def getfasta(bedfile,
             ref_genome,
             output_file,):
    """
    Extracts DNA sequences into a fasta file based on feature coordinates.
    """

    cmd_getfasta = "bedtools getfasta " \
                   "-fi %s.fa " \
                   "-bed %s " \
                   "-fo %s " \
                   "-tab" % (ref_genome,
                             bedfile,
                             output_file)
    return cmd_getfasta

