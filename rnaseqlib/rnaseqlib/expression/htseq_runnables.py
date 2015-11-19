"""
Collection of functions to execute htseq commands
"""


def htseq_count(uniq_sample_id,
                output_dir,
                outfile,
                bamfile,
                stranded,
                order,
                minaqual,
                mode,
                gtf):
    """
    Count reads in features with htseq-count
    :param uniq_patient_id: name says it all
    :param bamfile: bam input file
    :param stranded: Is the library stranded ? [yes,no]
    :param order:
    :param minaqual: min quality to filter for [default:10]
    :param mode: [union, intersection-strict, intersection-nonempty]
    :param gtf: gtf file
    :return:command to execute
    """

    cmd = "htseq-count"
    cmd += " --format bam"  # input is bamformat
    cmd += " --order {0}".format(order)
    cmd += " --stranded {0}".format(stranded)
    cmd += " -a {0}".format(minaqual)
    cmd += " --mode {0}".format(mode)
    #cmd += " --samout {0}/{1}".format(output_dir,uniq_sample_id) + ".sam"
    cmd += " {0} {1}".format(bamfile, gtf)
    cmd += " >{0}".format(outfile)

    return cmd