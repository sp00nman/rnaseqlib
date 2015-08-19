"""
Collection of functions that execute external software or UNIX commands.
"""


def varcall_bamfo(input_file,
                  output_file_gatk,
                  ref_genome):
    """
    Variant calling with bamfo (author: Tomas Konopka).
    """
    cmd_bamfo = "bamfo callvariants " + " \\\n" \
                + "--bam " + input_file + " \\\n" \
                + "--output " + output_file_gatk + " \\\n" \
                + "--genome " + ref_genome

    return cmd_bamfo
