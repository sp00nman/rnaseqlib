

def run_cmd(msg, cmd):
    logging.info(msg)
    logging.debug(cmd)
    status = 0
    if not args.debug:
        status = system(cmd)
        if status != 0:
            logging.warning("command '%s' returned non-zero "
                            "status: %d'" % (cmd, status))
    return status


def extract(input_file, sample_dir, project_name, region, output_dir):

    input_file = sample_dir + "/" + input_file
    output_file = output_dir + "/" + project_name + "_extract.bam"
    msg_extract = "Extract region: " + region
    cmd_extract = "samtools view -b -h %s %s >%s " % (input_file,
                                                      region,
                                                      output_file)
    return msg_extract, cmd_extract


def reorder(project_name, output_dir, ref_genome):

    input_file = output_dir + "/" + project_name + "_extract.bam"
    output_file = output_dir + "/" + project_name + "_reorder.bam"
    msg_replace = "Reorder readgroups."
    cmd_replace = "java -jar $NGS_PICARD/ReorderSam.jar " \
                  "INPUT=%s "\
                  "OUTPUT=%s " \
                  "REFERENCE=%s" % (input_file, output_file, ref_genome)
    return msg_replace, cmd_replace


def sort_bam(project_name, output_dir):
    """
    Sort sam ? bam file by coordinate.
    :param project_name: name of project (given by user)
    :param output_dir: where the output files should be written
    :return: message to be logged & command to be executed; type str
    """

    input_file =  output_dir + "/" + project_name + "_reorder.bam"
    output_file = output_dir + "/" + project_name + "_sorted.bam"
    msg_sort = "Sort bam file (by coordinate)."
    cmd_sort = "java -jar $NGS_PICARD/SortSam.jar " \
               "INPUT=%s " \
               "OUTPUT=%s " \
               "SORT_ORDER=coordinate" % (input_file, output_file)
    return msg_sort, cmd_sort


def replace_readgroups(project_name, output_dir):

    input_file = output_dir + "/" + project_name + "_sorted.bam"
    output_file = output_dir + "/" + project_name + "_replace.bam"
    msg_replace = "Reorder readgroups. "
    cmd_replace = "java -jar $NGS_PICARD/AddOrReplaceReadGroups.jar " \
                  "I=%s " \
                  "O=%s " \
                  "SO=coordinate " \
                  "RGID=1 " \
                  "RGLB=Lib1 " \
                  "RGPL=illumina " \
                  "RGPU=hiseq2000 " \
                  "RGSM=%s" % (input_file, output_file, project_name)
    return msg_replace, cmd_replace


def remove_duplicates(project_name, output_dir):
    """
    Remove duplicate reads.
    :param project_name: name of project (given by user)
    :param output_dir: where the output files should be written
    :return: message to be logged & command to be executed; type str
    """

    input_file = output_dir + "/" + project_name + "_replace.bam"
    output_file = output_dir + "/" + project_name + "_markduplicates.bam"
    output_metrics = output_dir + "/" + project_name + "_duplicates.metrics.txt"
    msg_rmdup = "Remove duplicate reads. "
    cmd_rmdup = "java -jar $NGS_PICARD/MarkDuplicates.jar " \
                "INPUT=%s " \
                "OUTPUT=%s " \
                "METRICS_FILE=%s " \
                "REMOVE_DUPLICATES=true" % (input_file, output_file,
                                            output_metrics)
    return msg_rmdup, cmd_rmdup


def index_bam(project_name, output_dir):
    """
    Index bam alignment file.
    :param project_name: name of project (given by user)
    """

    input_file = output_dir + "/" + project_name + "_markduplicates.bam"
    output_file = output_dir + "/" + project_name + "_markduplicates.bai"
    msg_indexbam = "Index bam file with samtools."
    cmd_indexbam = "samtools index %s %s" % (input_file, output_file)
    
    return msg_indexbam, cmd_indexbam


def splitntrim(project_name, output_dir, ref_genome):
    """
    SplitNCigarReads developed specially for RNAseq,
    which splits reads into exon segments (getting
    rid of Ns but maintaining grouping information)
    and hard-clip any sequences overhanging into the
    intronic regions.
    :param project_name: name of project (given by user
    :param output_dir: where the output files should be written
    :param ref_genome: reference genome  (.fa)
    :return: message to be logged & command to be executed; type str
    """

    input_file = output_dir + "/" + project_name + "_markduplicates.bam"
    output_file = output_dir + "/" + project_name + "_splitntrim.bam"
    msg_splitntrim = "Splitntrim. "
    cmd_splitntrim = "java -Xmx6g -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                     "-T SplitNCigarReads " \
                     "-R %s " \
                     "-I %s " \
                     "-o %s " \
                     "-rf ReassignOneMappingQuality " \
                     "-RMQF 255 " \
                     "-RMQT 60 " \
                     "-U ALLOW_N_CIGAR_READS" % (ref_genome, input_file,
                                                 output_file)
    return msg_splitntrim, cmd_splitntrim


def bqsr(project_name, output_dir, ref_genome):
    """
    We do recommend running base recalibration (BQSR). Even though the 
    effect is also marginal when applied to good quality data, it can 
    absolutely save your butt in cases where the qualities have 
    systematic error modes.
    Both steps 4 and 5 are run as described for DNAseq (with the same 
    known sites resource files), without any special arguments. Finally, 
    please note that you should NOT run ReduceReads on your RNAseq data. 
    The ReduceReads tool will no longer be available in GATK 3.0.
    """

    input_file = output_dir + "/" + project_name + "_splitntrim.bam"
    output_file = output_dir + "/" + project_name + "_bqsr.bam"
    msg_splitntrim = "Base recalibration. "
    cmd_splitntrim = "java -Xmx6g -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                     "-T PrintReads " \
                     "-R %s " \
                     "-I %s " \
                     "-BSQR recalibration_report.grp " \
                     "-o %s " % (ref_genome, input_file, output_file)
    
    return (msg_splitntrim, cmd_splitntrim)


def varcall_bamfo(project_name, output_dir, ref_genome):

    input_file = output_dir + "/" + project_name + "_bqsr.bam"
    output_file_gatk = output_dir + "/" + project_name + "_varcall_bamfo.vcf"
    msg_bamfo = "Bamfo variant calling."
    cmd_bamfo = "bamfo callvariants " + " \\\n" \
                + "--bam " + input_file + " \\\n" \
                + "--output " + output_file_gatk + " \\\n" \
                + "--genome " + ref_genome
    return msg_bamfo, cmd_bamfo


def varcall_samtools(project_name, output_dir, ref_genome):

    input_file = output_dir + "/" + project_name + "_bqsr.bam"
    output_file_samtools = output_dir + "/" + project_name + "_varcall_samtools.vcf"
    msg_samtools = "Samtools variant calling."
    cmd_samtools = "samtools mpileup " \
                   "-C 50 " \
                   "-uf %s " \
                   "%s " \
                   "| bcftools view -vcg - > %s" % (ref_genome, input_file,
                                                    output_file_samtools)
    return msg_samtools, cmd_samtools


def varcall_gatk(project_name, output_dir, ref_genome):
    """
    http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq
    'We have added some functionality to the variant calling code which
    will intelligently take into account the information about intron-exon
    split regions that is embedded in the BAM file by SplitNCigarReads.
    In brief, the new code will perform dangling head merging operations
    and avoid using soft-clipped bases (this is a temporary solution) as
    necessary to minimize false positive and false negative calls. To invoke
    this new functionality, just add -dontUseSoftClippedBases to your regular
    HC command line. Note that the -recoverDanglingHeads argument which was
    previously required is no longer necessary as that behavior is now enabled
    by default in HaplotypeCaller. Also, we found that we get better results
    if we lower the minimum phred-scaled confidence threshold for calling
    variants on RNAseq data, so we use a default of 20 (instead of 30 in
    DNA-seq data).'
    :param project_name:
    :param output_dir:
    :param ref_genome:
    :return:
    """
    input_file = output_dir + "/" + project_name + "_bqsr.bam"
    output_file_gatk = output_dir + "/" + project_name + "_varcall_gatk.vcf"
    msg_varcall_gatk = "Call variants (gatk)."
    cmd_varcall_gatk = "java -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                       "-T HaplotypeCaller " \
                       "-R %s " \
                       "-I %s " \
                       "-dontUseSoftClippedBases " \
                       "-stand_call_conf 20.0 " \
                       "-stand_emit_conf 20.0 " \
                       "-o %s" % (ref_genome, input_file, output_file_gatk)

    return msg_varcall_gatk, cmd_varcall_gatk


def variant_filtering(project_name, output_dir, ref_genome):
    """
    To filter the resulting callset, you will need to apply hard filters,
    as we do not yet have the RNAseq training/truth resources that would
    be needed to run variant recalibration (VQSR).
    We recommend that you filter clusters of at least 3 SNPs that are within
    a window of 35 bases between them by adding -window 35 -cluster 3 to your
    command. This filter recommendation is specific for RNA-seq data.
    As in DNA-seq, we recommend filtering based on Fisher Strand values
    (FS > 30.0) and Qual By Depth values (QD < 2.0).
    :param project_name:
    :param output_dir:
    :param ref_genome:
    :return:
    """

    input_file = output_dir + "/" + project_name + "_varcall_gatk.vcf"
    output_file = output_dir + "/" + project_name + "_filtering_gatk.bam"
    msg_filter = "Filtering."
    cmd_filter = "java -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                 "-T VariantFiltration " \
                 "-R %s " \
                 "-V %s " \
                 "-window 35 " \
                 "-cluster 3 " \
                 "-filterName FS " \
                 "-filter \"FS > 30.0\" " \
                 "-filterName QD " \
                 "-filter \"QD < 2.0\" " \
                 "-o %s" % (ref_genome, input_file, output_file)

    return msg_filter, cmd_filter