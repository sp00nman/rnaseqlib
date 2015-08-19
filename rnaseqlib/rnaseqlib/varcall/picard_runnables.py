"""
Collection of functions that execute external software or UNIX commands.
"""


def reorder_sam(inbamfile,
                outbamfile,
                genome_path,
                heap_mem):
    """
    Reorder BAM file.
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file
    :param genomes: path to genome
    :return: Command to be executed; type str
    """

    cmd_reorder = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                  "-jar $NGS_PICARD/ReorderSam.jar " \
                  "INPUT=%s " \
                  "OUTPUT=%s " \
                  "REFERENCE=%s.fa" % (heap_mem,
                                       inbamfile,
                                       outbamfile,
                                       genome_path)
    return cmd_reorder


def sort_bam(inbamfile,
             outbamfile,
             heap_mem,
             sort_order="coordinate"):
    """
    Sort BAM file by variable
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file (BAM formatted and sorted)
    :param sort_order: sort by (default: coordinate)
    :return: Command to be executed; type str
    """

    cmd_sort = "java -%s -Djava.io.tmpdir=$TMPDIR " \
               "-jar $NGS_PICARD/SortSam.jar " \
               "INPUT=%s " \
               "OUTPUT=%s " \
               "SORT_ORDER=%s" % (heap_mem,
                                  inbamfile,
                                  outbamfile,
                                  sort_order)
    return cmd_sort


def replace_readgroups(input_file,
                       output_file,
                       project_name,
                       heap_mem):
    """
    Replace readgroups (with some dummy variables.)
    """

    cmd_replace = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                  "-jar $NGS_PICARD/AddOrReplaceReadGroups.jar " \
                  "I=%s " \
                  "O=%s " \
                  "SO=coordinate " \
                  "RGID=1 " \
                  "RGLB=Lib1 " \
                  "RGPL=illumina " \
                  "RGPU=hiseq2000 " \
                  "RGSM=%s" % (heap_mem,
                               input_file,
                               output_file,
                               project_name)
    return cmd_replace


def remove_duplicates(inbamfile,
                      outbamfile,
                      metrics_file,
                      heap_mem):
    """
    Remove duplicate reads.
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file
    :param metrics_file: file with summary statistics about
    :return: Command to be executed; type str
    """

    cmd_rmdup = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                "-jar $NGS_PICARD/MarkDuplicates.jar " \
                "INPUT=%s " \
                "OUTPUT=%s " \
                "METRICS_FILE=%s " \
                "REMOVE_DUPLICATES=true" % (heap_mem,
                                            inbamfile,
                                            outbamfile,
                                            metrics_file)
    return cmd_rmdup
