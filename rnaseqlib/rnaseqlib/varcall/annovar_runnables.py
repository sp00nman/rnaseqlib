"""
Collection of functions that execute external software or UNIX commands.
"""


def run_annovar(vcf_file,
                protocol,
                output_file,
                operation,
                buildversion,
                annovar_dir,):
    """
    Annotate VCF file with annovar_150322.
    """

    cmd_annovar = "table_annovar.pl " \
                  "--protocol %s " \
                  "--operation %s " \
                  "--outfile %s " \
                  "--buildver %s " \
                  "--remove " \
                  "--otherinfo " \
                  "--nastring . " \
                  "--vcfinput " \
                  "%s " \
                  "%s" % (protocol,
                          operation,
                          output_file,
                          buildversion,
                          vcf_file,
                          annovar_dir)

    return cmd_annovar


def convert2annovar(input_file,
                    output_file):
    """
    Annotate vcf file with annovar.
    """

    cmd_annovar = "convert2annovar.pl " \
                  "--format vcf4 " \
                  "--includeinfo %s >" \
                  "%s.annovar" % (input_file,
                                 output_file)
    return cmd_annovar


def dbsnp_filter(dbtype,
                 buildversion,
                 input_file,
                 annovar_dir):
    """
    Filter for dbsnp database.
    """

    cmd_dbsnp_filter = "annotate_variation.pl " \
                        "--filter "\
                        "--dbtype %s " \
                        "--buildver %s " \
                        "%s %s" % (dbtype,
                                   buildversion,
                                   input_file,
                                   annovar_dir)
    return cmd_dbsnp_filter


def gene_annotation(buildversion,
                    input_file,
                    annovar_dir):
    """
    Annotate variants.
    """

    cmd_dbsnp_filter = "annotate_variation.pl " \
                        "--buildver %s " \
                        "%s %s" % (buildversion,
                                   input_file,
                                   annovar_dir)
    return cmd_dbsnp_filter

