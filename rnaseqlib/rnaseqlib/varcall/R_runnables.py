"""
Collection of functions to execute R script files.
"""


def plot_patient_mutation_matrix(
        clinical_info,
        mutation_info,
        output_file,
        dn
):
    """
    Plot patient mutation matrix
    :param clinical_info: clinical information from patients (eg. diagnosis,
    mutation status, gender, etc...
    :param mutation_info: mutations per patients and other metrics like
    allele frequency, etc...
    :param dn: execution directory
    :return:
    """
    cmd_plot = "Rscript --vanilla " \
               + dn + "/" \
               + "plot_patient_mutation_matrix.R %s %s %s" % (mutation_info,
                                                              clinical_info,
                                                              output_file)
    return cmd_plot