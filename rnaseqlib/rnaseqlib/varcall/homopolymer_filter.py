"""
Filter for variants that fall within homopolymer sites.
"""

import vcf
import re
from rnaseqlib.utils import tools as ts


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
            # zb AAAA [A] AAAA ; [A]=variant
            # 4 nucleotides left and right from variant position
            out_file.writelines(
                record.CHROM + "\t"
                + str(record.POS - (nucleo_num +1)) + "\t"
                + str(record.POS + nucleo_num) + "\n"
            )

    finally:
            vcf_file.close()
            out_file.close()

    return None


def annotate_variants(coord_seq,
                      input_file,
                      output_file,
                      nucleo_num):

    coord_dic = ts.load_dictionary(
        input_file=coord_seq,
        sep='\t'
    )

    in_handle = open(input_file, 'r')
    out_handle = open(output_file, 'w')

    #TODO: print header first

    try:
        vcf_reader = vcf.Reader(in_handle)
        vcf_writer = vcf.Writer(out_handle, vcf_reader)

        for record in vcf_reader:
            key = record.CHROM + ":" \
                  + str(record.POS - (nucleo_num+1)) + "-" \
                  + str(record.POS + nucleo_num)

            seq = coord_dic[key]
            seq_homopolymer = seq[nucleo_num] * nucleo_num

            if (re.match(seq_homopolymer,
                         seq[0:nucleo_num-1]) or \
                        re.match(seq_homopolymer,
                                 seq[nucleo_num+1:len(seq)]
                        )
            ):
                record.FILTER.append('HRun')

            vcf_writer.write_record(record)

    finally:
            in_handle.close()
            out_handle.close()
