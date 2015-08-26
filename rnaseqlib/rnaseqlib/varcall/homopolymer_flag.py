"""
Flag variants that fall within homopolymer sites (HRun) or
near homopolymer sites (nHRun).
"""

import vcf
import re
from rnaseqlib.utils import tools as ts


def extract_pos2bed(input_file,
                    bedfile,
                    nucleo_num,
                    dis):
    """
    Extract chromosome start end from vcf and write as BED file. Add a
    sepcified number of nucleotides left and right to the position.
    :param input_file: vcf file
    :param bedfile: Write output to BED file
    :param nucleo_num: Number of nucleotide left and right of position
    :param dis: Repeat search dis nucleotides away from variant
    :return: None
    """

    vcf_file = open(input_file, 'r')
    out_file = open(bedfile, 'w')
    # extracts +1 nucleotide from left and right for nHRun filter
    nhrun = nucleo_num + dis

    try:
        vcf_reader = vcf.Reader(vcf_file)

        for record in vcf_reader:
            # zb AAAA [A] AAAA ; [A]=variant
            # 4 nucleotides left and right from variant position
            out_file.writelines(
                record.CHROM + "\t"
                + str(record.POS - (nhrun +1)) + "\t"
                + str(record.POS + nhrun) + "\n"
            )

    finally:
            vcf_file.close()
            out_file.close()

    return None


def annotate_variants(coord_seq,
                      input_file,
                      output_file,
                      nucleo_num,
                      dis):

    nhrun = nucleo_num + dis

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
                  + str(record.POS - (nhrun+1)) + "-" \
                  + str(record.POS + nhrun)

            seq = coord_dic[key]

            for i in range(0,dis+1):
                # repeat search dis nucleotides away
                for lr in [dis,-dis]:
                    seq_homopolymer = seq[nucleo_num + i + lr] * nucleo_num
                    hit = match_homopolymer(seq_homopolymer,
                                            seq,
                                            nucleo_num,
                                            dis,
                                            i)
                    print key, seq, seq[nucleo_num + i + lr], seq_homopolymer, dis, i
                    if (hit):
                        if i==0:
                            record.FILTER.append('HRun')
                            # leave loop
                            break
                        else:
                            #TODO: would make sense to incorporate
                            # a [dis]nHRun filter
                            record.FILTER.append('nHRun')
                    # leave loop after first match was found
                    break

            vcf_writer.write_record(record)

    finally:
        in_handle.close()
        out_handle.close()


def match_homopolymer(seq_homopolymer,
                      seq,
                      nucleo_num,
                      dis,
                      i):
    return re.match(seq_homopolymer,
                    seq[0:(nucleo_num + dis +i) - 1]) or \
        re.match(seq_homopolymer,
                 seq[(nucleo_num + dis + i) + 1:len(seq)])
