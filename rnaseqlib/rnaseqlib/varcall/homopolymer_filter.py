"""
Filter for variants that fall within homopolymer sites.
"""

import vcf
import re
from rnaseqlib.utils import tools as ts
from rnaseqlib.varcall import runnables as rb


def extract_pos2bed(input_file,
                    bedfile,
                    nucleo_num):
    """
    Extract chromosome start end from vcf and write as BED file. Add a
    sepcified number of nucleotides left and right to the position.
    :param input_file: vcf file
    :param bedfile: Write output to BED file
    :param nucleo_num: Number of nucleotide left and right of position
    :return: None
    """

    vcf_file = open(input_file, 'r')
    out_file = open(bedfile, 'w')

    try:
        vcf_reader = vcf.Reader(vcf_file)

        for record in vcf_reader:
            out_file.writelines(
                record.CHROM + "\t"
                + str(record.POS - nucleo_num) + "\t"
                + str(record.POS + nucleo_num) + "\n"
            )

    finally:
            vcf_file.close()
            out_file.close()

    return None


def get_coordinates(bedfile,
                    ref_genome,
                    output_file):

    cmd = rb.bedtools_getfasta(
        bedfile=bedfile,
        ref_genome=ref_genome,
        output_file=output_file
    )

    return cmd


def annotate_variants(coordinates,
                      input_file,
                      output_file,
                      nucleo_num):

    coord_dic = ts.load_dictionary(
        input_file=coordinates,
        sep='\t'
    )

    vcf_file = open(input_file, 'r')
    vcf_writer = vcf.Writer(open(output_file, 'w'))

    #TODO: print header first

    try:
        vcf_reader = vcf.Reader(vcf_file)

        for record in vcf_reader:
            key = record.CHROM + ":" \
                  + str(record.POS - nucleo_num) + "-" \
                  + str(record.POS + nucleo_num)

            seq = coord_dic[key]
            seq_homopolymer = seq[nucleo_num+1]*nucleo_num
            if (re.match(seq_homopolymer, seq[0:nucleo_num]) or \
                        re.match(seq_homopolymer, seq[nucleo_num+1])):
                record.FILTER.append('HRun')
                vcf_writer.write_record(record)

    finally:
            vcf_file.close()

