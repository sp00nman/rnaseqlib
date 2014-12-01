#!/usr/bin/env python

# author fschischlik
# rna seq variant calling; for now limited to a specific region
# gatk OR bamfo OR samtools OR all

import argparse
import re
from os import system
import os
import logging


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


def extract(input_file, project_name, region, output_dir):

    output_file = output_dir + "/" + project_name + "_extract.bam"
    msg_extract = "Extract region: " + region
    cmd_extract = "samtools view -b -h %s %s >%s " % (input_file,
                                                      region,
                                                      output_file)
    return msg_extract, cmd_extract


def reorder(input_file, project_name, output_dir, ref_genome):

    output_file = output_dir + "/" + project_name + "_reorder.bam"
    msg_replace = "Reorder readgroups."
    cmd_replace = "java -jar $NGS_PICARD/ReorderSam.jar " \
                  "INPUT=%s "\
                  "OUTPUT=%s " \
                  "REFERENCE=%s" % (input_file, output_file, ref_genome)
    return msg_replace, cmd_replace


def sort_bam(input_file, project_name, output_dir):
    """
    Sort sam ? bam file by coordinate.
    :param project_name: name of project (given by user)
    :param output_dir: where the output files should be written
    :return: message to be logged & command to be executed; type str
    """

    output_file = output_dir + "/" + project_name + "_sorted.bam"
    msg_sort = "Sort bam file (by coordinate)."
    cmd_sort = "java -jar $NGS_PICARD/SortSam.jar " \
               "INPUT=%s " \
               "OUTPUT=%s " \
               "SORT_ORDER=coordinate" % (input_file, output_file)
    return msg_sort, cmd_sort


def replace_readgroups(input_file, project_name,
                       output_dir):

    output_file = output_dir + "/" + project_name + "_replace.bam"
    msg_replace = "Reorder readgroups. "
    cmd_replace = "java -jar $NGS_PICARD/AddOrReplaceReadGroups.jar " \
                  "I=%s " \
                  "O=%s " \
                  "SO=coordinate" \
                  "RGID=1" \
                  "RGLB=Lib1" \
                  "RGPL=illumina" \
                  "RGPU=hiseq2000" \
                  "RGSM=%s" % (input_file, output_file, project_name)
    return msg_replace, cmd_replace


def remove_duplicates(input_file, project_name,
                      output_dir):
    """
    Remove duplicate reads.
    :param project_name: name of project (given by user)
    :param output_dir: where the output files should be written
    :return: message to be logged & command to be executed; type str
    """

    output_file = output_dir + "/" + project_name + "_markduplicates.bam"
    msg_rmdup = "Remove duplicate reads. "
    cmd_rmdup = "java -jar $NGS_PICARD/MarkDuplicates.jar " \
                "INPUT=%s " \
                "OUTPUT=%s " \
                "METRICS_FILE=%s.duplicates.metrics.txt " \
                "REMOVE_DUPLICATES=true" % (input_file, output_file,
                                            project_name)
    return msg_rmdup, cmd_rmdup


def splitntrim(input_file, project_name, output_dir, ref_genome):

    output_file = output_dir + "/" + project_name + "_split.bam"
    msg_splitntrim = "Splitntrim. "
    cmd_splitntrim = "java -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                     "-T SplitNCigarReads" \
                     "-R %s" \
                     "-I %s" \
                     "-o %s" \
                     "-rf ReassignOneMappingQuality" \
                     "-RMQF 255" \
                     "-RMQT 60" \
                     "-U ALLOW_N_CIGAR_READS" % (ref_genome, input_file,
                                                 output_file)
    return msg_splitntrim, cmd_splitntrim


def variant_calling(input_file, project_name, output_dir, ref_genome):

    #bamfo

    #samtools

    #gatk
    output_file_gatk = output_dir + "/" + project_name + "_varcall_gatk.vcf"
    msg_varcall_gatk = "Call variants (gatk)."
    cmd_varcall_gatk = "java -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                       "-T HaplotypeCaller " \
                       "-R %s " \
                       "-I %s " \
                       "-dontUseSoftClippedBases " \
                       "-stand_call_conf 20.0" \
                       "-stand_emit_conf 20.0" \
                       "-o %s" % (ref_genome, input_file, output_file_gatk)

    return msg_varcall_gatk, cmd_varcall_gatk


def variant_filtering(input_file, project_name, output_dir, ref_genome):

    output_file = output_dir + "/" + project_name + "_filtering.bam"
    msg_filter = "Filtering."
    cmd_filter = "java -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                 "-T VariantFiltration " \
                 "-R %s " \
                 "-V %s " \
                 "-window 35 " \
                 "-cluster 3 " \
                 "-filterName FS " \
                 "-filter \"FS > 30.0\"" \
                 "-filterName QD " \
                 "-filter \"QD < 2.0\"" \
                 "-o %s" % (ref_genome, input_file, output_file)

    return msg_filter, cmd_filter


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genetic screen workflow '
                                                 '0.0.1')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False, default="all",
                        choices=["all", "extract", "reorder", "duplicates",
                                 "splitntrim", "realignment", "recalibration",
                                 "varcall", "filter"],
                        help='Limit job submission to a particular '
                             'Analysis stage.')
    parser.add_argument('--project_name', required=False, type=str,
                        help="name of the project")
    parser.add_argument('--input_file', required=False, type=str,
                        help="name of the sample file")
    parser.add_argument('--sample_dir', required=False, type=str,
                        help="sample directory")
    parser.add_argument('--exec_dir', required=False, type=str,
                        help="exec_dir")
    parser.add_argument('--output_dir', required=False, type=str,
                        help="output_dir")
    parser.add_argument('--ref_genome', required=False, type=str,
                        help="reference genome")
    parser.add_argument('--region', required=False, type=str,
                        help="region eg. 1:23842983-23843983")

    # parse arguments, set defaults
    args = parser.parse_args()
    home_dir = os.getenv("HOME")
    if not args.output_dir:
        args.output_dir = os.getcwd()
    if not args.sample_dir:
        args.sample_dir = os.getcwd()
    if not args.exec_dir:
        args.exec_dir = home_dir + "/src"
    if not args.ref_genome:
        args.defuse_ref = home_dir + "/ref_genome"

    # logging
    # create log file
    logfile_name = args.project_name + ".log"
    logging.basicConfig(filename=logfile_name,
                        format='%(levelname)s: %(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)

    # start analysis workflow & logging
    logging.info("RNAseq variant calling (region specific)")

    # start workflow
    if re.search(r"all|extract", args.stage):
        (msg, cmd) = extract(args.input_file, args.project_name,
                             args.region, args.output_dir)
        status = run_cmd(msg, cmd)

    if re.search(r"all|reorder", args.stage):
        (msg, cmd) = reorder(args.input_file, args.project_name,
                             args.output_dir, args.ref_genome)
        status = run_cmd(msg, cmd)

    if re.search(r"all|sort_bam", args.stage):
        (msg, cmd) = reorder(args.input_file, args.project_name,
                             args.output_dir)
        status = run_cmd(msg, cmd)

    if re.search(r"all|replace_readgroups", args.stage):
        (msg, cmd) = reorder(args.input_file, args.project_name,
                             args.output_dir)
        status = run_cmd(msg, cmd)

    if re.search(r"all|remove_duplicates", args.stage):
        (msg, cmd) = reorder(args.input_file, args.project_name,
                             args.output_dir)
        status = run_cmd(msg, cmd)
