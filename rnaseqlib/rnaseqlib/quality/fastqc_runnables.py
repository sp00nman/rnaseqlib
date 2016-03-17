"""
Collection of functions that execute external software or UNIX commands.
"""


def fastqc(input_file,
           output_dir):
    """
    Run fastqc on BAM file. Uses sequence file format detection.
    :param input_file: BAM file
    :param output_dir: Create all output files in the specified output directory
    :return: str object; Command to be executed
    """

    cmd_fastqc = "fastqc"
    cmd_fastqc += " input_file"
    cmd_fastqc += " -o " + output_dir

    return cmd_fastqc
