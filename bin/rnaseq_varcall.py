#!/usr/bin/env python

import argparse
import re
import os
import logging

from rnaseqlib.utils import tools as ts
from rnaseqlib.utils import parse_annovar as pa
from rnaseqlib.varcall import samtools_runnables as samtools_rb
from rnaseqlib.varcall import gatk_runnables as gatk_rb
from rnaseqlib.varcall import picard_runnables as picard_rb
from rnaseqlib.varcall import star_runnables as star_rb
from rnaseqlib.varcall import annovar_runnables as annovar_rb
from rnaseqlib.varcall import bedtools_runnables as bedtools_rb
from rnaseqlib.varcall import bamfo_runnables as bamfo_rb
from rnaseqlib.varcall import homopolymer_flag as hf
from rnaseqlib.varcall import near_indel_flag as nif
from rnaseqlib.varcall import filter_vcf as fv


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genetic screen workflow '
                                                 '0.0.1')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'analysis stage.'
                        '[all,alignment,star2pass,extract,duplicates,index,'
                        'splitntrim,indel,bqsr,bamfo,samtools,gatk,gatk_flag,'
                        'hrun_flag,nind_flag,annovar,selection]')
    parser.add_argument('--project_name', required=False, type=str,
                        help="name of the project")
    parser.add_argument('--read1', required=False, type=str,
                        help="For paired alignment, forward read.")
    parser.add_argument('--read2', required=False, type=str,
                        help="For paired alignment, reverse read.")

    parser.add_argument('--sample_dir', required=False, type=str,
                        help="Path to sample directory")
    parser.add_argument('--output_dir', required=False, type=str,
                        help="Path to output directory.")

    # aligner options
    parser.add_argument('--star_genome', required=False, type=str,
                        help="Genome directory of star aligner.")
    parser.add_argument('--ref_genome', required=False, type=str,
                        help="reference genome")
    parser.add_argument('--gtf', required=False, type=str,
                        help="gtf file annotation")

    # input files for indel realignment
    # for now at least 2 known indel files have to be provided
    parser.add_argument('--known_indels_1000g', required=False, type=str,
                        help="List of known indel sites [eg. 1000G_phase1.indels.b37.vcf]. ")
    parser.add_argument('--known_indels_mills', required=False, type=str,
                        help="List of known indel sites. [eg. Mills_and_1000G_gold_standard.indels.b37.vcf] ")

    # input files for base recalibration
    parser.add_argument('--dbsnp', required=False, type=str,
                        help="List of known snps.")

    # options for region specific variant calling
    parser.add_argument('--sample_file', required=False, type=str,
                        help="For region variant calling.")
    parser.add_argument('--region', required=False, type=str, default=False,
                        help="region eg. 20:30946147-31027122")

    # annovar specific options
    parser.add_argument('--annovar', required=False, type=str,
                        help="Annotate variant with annovar.")

    # hardware specific options
    parser.add_argument('--num_cpus', dest='num_cpus', required=False,
                        help='Number of cpus.')
    parser.add_argument('--heap_mem', dest='heap_mem', required=False,
                        help='Maximum heap size provided to Java. [Xmx[num]g]')


    # defaults
    args = parser.parse_args()
    home_dir = os.getenv("HOME")
    #TODO: implement this as optional
    read_length = 100

    if not args.output_dir:
        args.output_dir = os.getcwd()
    if not args.sample_dir:
        args.sample_dir = os.getcwd()
    if not args.ref_genome:
        args.defuse_ref = home_dir + "/ref_genome"
    if not args.star_genome:
        args.star_genome = home_dir + "/star_genome"
    if not args.annovar:
        args.annovar = home_dir + "/" + "src" + "/" + "annovar" + "/" + "humandb"
    if not args.num_cpus:
        args.num_cpus = "1"
    if not args.heap_mem:
        args.heap_mem = "Xmx6g"

    #

    # set project directory
    project_dir = args.output_dir + "/" + args.project_name

    # create directory structure
    ts.create_output_dir(args.output_dir, args.project_name)
    sub_dir = args.output_dir + "/" + args.project_name
    ts.create_output_dir(sub_dir, "star_2pass")

    # create log file
    logfile_name = args.output_dir + "/" + args.project_name + "/" \
                   + args.project_name + ".log"
    logging.basicConfig(filename=logfile_name,
                        format='%(levelname)s: %(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)

    # start analysis workflow & logging
    logging.info("Initiated RNA-seq variant calling workflow.")

    # load dictionary for file extensions
    data_dir = str(os.path.dirname(os.path.realpath(__file__)).strip('bin')) \
               + "data" + "/"
    file_ext = ts.load_dictionary(data_dir + 'file_extension.txt')
    # load dictionary for stdout messages
    stdout_msg = ts.load_dictionary(data_dir + 'stdout_message.txt')

    # start workflow
    if re.search(r"all|alignment", args.stage):

        cmd = star_rb.rnaseq_align(
            star_genome=args.star_genome,
            gtf=args.gtf,
            read1=args.read1,
            read2=args.read2,
            num_cpus=args.num_cpus,
            outfile_prefix=project_dir + "/"
                           + args.project_name + "."
        )

        status = ts.run_cmd(
            message=stdout_msg['alignment'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|star2pass", args.stage):
        cmd = star_rb.star_index(
            novel_ref=project_dir + "/"
                      + "star_2pass",
            genome=args.ref_genome,
            firstroundalignment=project_dir + "/"
                                + args.project_name + "."
                                + file_ext['sj_out_tab'],
            sjdb_overhang=read_length-1,
            num_cpus=args.num_cpus,
            outfile_prefix=project_dir + "/"
                           + args.project_name + "."
        )

        status = ts.run_cmd(
            message=stdout_msg['alignment_index'],
            command=cmd,
            debug=args.debug
        )

        cmd = star_rb.rnaseq_align(
            star_genome=project_dir + "/"
                       + "star_2pass",
            gtf=args.gtf,
            read1=args.read1,
            read2=args.read2,
            num_cpus=args.num_cpus,
            outfile_prefix=project_dir + "/"
                          + args.project_name + "."
                          + "2pass_"
        )

        status = ts.run_cmd(
            message=stdout_msg['realign'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|extract", args.stage):

        if args.region is False:
            sample_file = project_dir + "/" \
                          + args.project_name + "." \
                          + file_ext['alignment_index']
        else:
            sample_file = args.sample_file

        cmd = samtools_rb.extract(
            input_file=sample_file,
            region=args.region,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['extract']
        )

        status = ts.run_cmd(
            message=stdout_msg['extract'],
            command=cmd,
            debug=args.debug
        )

        # reorder
        if args.region is False:
            sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['extract']
        else:
            if args.star2pass:
                sample_file = project_dir + "/" \
                              + args.project_name + "." \
                              + file_ext['star2pass']
            else:
                sample_file = project_dir + "/" \
                              + args.project_name + "." \
                              + file_ext['star_alignment']

        cmd = picard_rb.reorder_sam(
            inbamfile=sample_file,
            outbamfile=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['reorder'],
            genome_path=args.ref_genome)

        status = ts.run_cmd(
            message=stdout_msg['reorder'],
            command=cmd,
            debug=args.debug
        )

        # sort bam
        sample_file = project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['reorder']

        cmd = picard_rb.sort_bam(
            inbamfile=sample_file,
            outbamfile=project_dir + "/"
                        + args.project_name + "." \
                        + file_ext['sort'],
            sort_order="coordinate",
            heap_mem=args.heap_mem)

        status = ts.run_cmd(
            message=stdout_msg['sort'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|replace_rg", args.stage):

        if args.region is False:
            if re.search(r"star2pass", args.stage):
                sample_file = project_dir + "/" \
                            + args.project_name + "." \
                            + file_ext['star2pass']
            else:
                sample_file = project_dir + "/" \
                            + args.project_name + "." \
                            + file_ext['star_alignment']
        else:
            sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['sort']

        cmd = picard_rb.replace_readgroups(
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['replace_rg'],
            project_name=args.project_name,
            heap_mem=args.heap_mem)

        status = ts.run_cmd(
            message=stdout_msg['replace_rg'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|duplicates", args.stage):

        sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['replace_rg']

        cmd = picard_rb.remove_duplicates(
            inbamfile=sample_file,
            outbamfile=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['duplicates'],
            metrics_file=".duplicate_metrics.txt",
            heap_mem=args.heap_mem)

        status = ts.run_cmd(
            message=stdout_msg['duplicates'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|index", args.stage):

        sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['duplicates']

        cmd = samtools_rb.index_bam(
            inbamfile=sample_file)

        status = ts.run_cmd(
            message=stdout_msg['index'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|splitntrim", args.stage):

        sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['duplicates']

        cmd = gatk_rb.splitntrim(
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['splitntrim'],
            ref_genome=args.ref_genome,
            heap_mem=args.heap_mem)

        status = ts.run_cmd(
            message=stdout_msg['splitntrim'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|indel", args.stage):

        sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['splitntrim']

        cmd = gatk_rb.target_creator(
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['target'],
            ref_genome=args.ref_genome,
            known_1000g=args.known_indels_1000g,
            known_mills=args.known_indels_mills,
            heap_mem=args.heap_mem

        )

        status = ts.run_cmd(
            message=stdout_msg['target'],
            command=cmd,
            debug=args.debug
        )

        cmd = gatk_rb.indel_realignment(
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['indel'],
            ref_genome=args.ref_genome,
            target=project_dir + "/"
                   + args.project_name + "."
                   + file_ext['target'],
            known_1000g=args.known_indels_1000g,
            known_mills=args.known_indels_mills,
            heap_mem=args.heap_mem
        )

        status =ts.run_cmd(
            message=stdout_msg['indel'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|bqsr", args.stage):

        sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['indel']

        cmd = gatk_rb.analyse_covariation_patterns(
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['bqsrI'],
            ref_genome=args.ref_genome,
            known_1000g=args.known_indels_1000g,
            known_mills=args.known_indels_mills,
            dbsnp=args.dbsnp,
            heap_mem=args.heap_mem
        )

        status = ts.run_cmd(
            message=stdout_msg['bqsrI'],
            command=cmd,
            debug=args.debug
        )

        cmd = gatk_rb.analyse_covariation_patterns_2ndpass(
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['bqsrII'],
            ref_genome=args.ref_genome,
            known_1000g=args.known_indels_1000g,
            known_mills=args.known_indels_mills,
            dbsnp=args.dbsnp,
            heap_mem=args.heap_mem,
            recal=project_dir + "/"
                  + args.project_name + "."
                  + file_ext['bqsrI']
        )

        status = ts.run_cmd(
            message=stdout_msg['bqsrII'],
            command=cmd,
            debug=args.debug
        )

        cmd = gatk_rb.plot_recalibration(
            ref_genome=args.ref_genome,
            before=project_dir + "/"
                   + args.project_name + "."
                   + file_ext['bqsrI'],
            after=project_dir + "/"
                  + args.project_name + "."
                  + file_ext['bqsrII'],
            plot_name=project_dir + "/"
                     + args.project_name + "."
                     + file_ext['recal_plot'],
            heap_mem=args.heap_mem
        )

        status = ts.run_cmd(
            message=stdout_msg['plot_recal'],
            command=cmd,
            debug=args.debug
        )

        cmd = gatk_rb.bqsr(
            ref_genome=args.ref_genome,
            input_file=project_dir + "/" \
                       + args.project_name + "." \
                       + file_ext['indel'],
            bqsr=project_dir + "/"
                 + args.project_name + "."
                 + file_ext['bqsrI'],
            output_file=project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['bqsr'],
            heap_mem=args.heap_mem
        )

        status = ts.run_cmd(
            message=stdout_msg['bqsr'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|bamfo", args.stage):

        sample_file = project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['bqsr']

        cmd = bamfo_rb.varcall_bamfo(
            input_file=sample_file,
            output_file_gatk=project_dir + "/"
                            + args.project_name + "."
                            + file_ext['bamfo'],
            ref_genome=args.ref_genome
        )

        status = ts.run_cmd(
            message=stdout_msg['bamfo'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|samtools", args.stage):

        sample_file = project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['bqsr']

        cmd = samtools_rb.varcall_samtools(
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['samtools'],
            ref_genome=args.ref_genome)

        status = ts.run_cmd(
            message=stdout_msg['samtools'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|gatk", args.stage):

        sample_file = project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['bqsr']

        cmd = gatk_rb.varcall_gatk(
            input_file=sample_file,
            output_file_gatk=project_dir + "/"
                            + args.project_name + "."
                            + file_ext['gatk'],
            ref_genome=args.ref_genome,
            heap_mem=args.heap_mem
        )

        status = ts.run_cmd(
            message=stdout_msg['gatk'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|gatk_flag", args.stage):

        #only gatk
        sample_file = project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['gatk']

        cmd = gatk_rb.variant_filtering(
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['filtering'],
            ref_genome=args.ref_genome,
            heap_mem=args.heap_mem
        )

        status = ts.run_cmd(
            message=stdout_msg['filtering'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|hrun_flag", args.stage):

        sample_file = project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['filtering']

        hf.extract_pos2bed(
            input_file=sample_file,
            bedfile=project_dir + "/"
                    + args.project_name + "."
                    + file_ext['coordinates'],
            nucleo_num=4,
            dis=1
        )

        # get coordinates from vcf file
        cmd = bedtools_rb.getfasta(
            bedfile=project_dir + "/"
                    + args.project_name + "."
                    + file_ext['coordinates'],
            ref_genome=args.ref_genome,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['seq']
        )

        status = ts.run_cmd(
            message=stdout_msg['hrun'],
            command=cmd,
            debug=args.debug
        )

        hf.annotate_variants(
            coord_seq=project_dir + "/"
                      + args.project_name + "."
                      + file_ext['seq'],
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['hrun'],
            nucleo_num=4,
            dis=1
        )

    if re.search(r"all|nind_flag", args.stage):
        # near indel, indel matches with nindel_flag -->so stupid

        sample_file = project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['hrun']

        nif.flag_variants(
            input_file=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['nindel'],
            dis=5
        )

    if re.search(r"all|annovar", args.stage):

        sample_file = project_dir + "/" \
                     + args.project_name + "." \
                     + file_ext['nindel']

        cmd = annovar_rb.run_annovar(
            vcf_file=sample_file,
            protocol="ensGene,"
                     + "cytoBand,"
                     + "genomicSuperDups,"
                     + "snp129,"
                     + "snp142Common,"
                     + "1000g2014oct_all,"
                     + "1000g2014oct_afr,"
                     + "1000g2014oct_amr,"
                     + "1000g2014oct_eas,"
                     + "1000g2014oct_eur,"
                     + "1000g2014oct_sas,"
                     + "esp6500siv2_all,"
                     + "cosmic70,"
                     + "clinvar_20150330,"
                     + "ljb26_all,"
                     + "caddgt10,"
                     + "caddindel",
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['annotation'],
            operation="g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f",
            buildversion="hg19",
            annovar_dir=args.annovar
        )

        status = ts.run_cmd(
            message=stdout_msg['annotation'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|selection", args.stage):

        sample_file = project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['annovar']

        fv.filter_vcf(
            input_file=sample_file,
            file_prefix=project_dir + "/" \
                        + args.project_name,
            file_suffix=file_ext['inhouse']
        )
