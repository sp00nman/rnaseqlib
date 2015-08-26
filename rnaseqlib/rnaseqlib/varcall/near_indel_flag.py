"""
Flag variants that fall within 5 bp of an indel.
VCF file has to be sorted.
"""

import vcf


def flag_variants(input_file,
                  output_file,
                  dis):

    in_reader = open(input_file, 'r')
    in_walker = open(input_file, 'r')
    in_counter = open(input_file, 'r')
    out_handle = open(output_file, 'w')

    try:
        vcf_reader = vcf.Reader(in_reader)
        vcf_walker = vcf.Reader(in_walker)
        vcf_counter = vcf.Reader(in_counter)
        vcf_writer = vcf.Writer(out_handle, vcf_reader)
        store_records = [None]*2

        # count the number of entries
        # until I find a better solution
        count_entries = 0
        for record in vcf_counter:
            count_entries +=1

        # initiate the walker
        store_records[0] = vcf_walker.next()

        # set counter to 0
        counter = 0

        for record in vcf_reader:
            counter += 1

            if counter == 1:
                # process first record
                store_records[1] = vcf_walker.next()
                        #     if record.is_snp:
                if record.CHROM == store_records[1].CHROM and \
                            (store_records[1].POS - record.POS) <= dis and \
                            store_records[1].is_indel:
                        record.FILTER.append('nIndel')

            if counter == count_entries:
                # process last record
                if record.is_snp:
                    if record.CHROM == store_records[0].CHROM and \
                            (store_records[0].POS - record.POS) <= dis and \
                            store_records[0].is_indel:
                        record.FILTER.append('nIndel')

            elif 0 > counter < count_entries:
                # process everything else
                store_records[1] = vcf_walker.next()
                if record.is_snp:
                    if record.CHROM == store_records[0].CHROM and \
                            store_records[0].is_indel and \
                            (record.POS - store_records[0].POS) <= dis or \
                            record.CHROM == store_records[1].CHROM and \
                            store_records[1].is_indel and \
                            (store_records[1].POS - record.POS) <= dis:
                         record.FILTER.append('nIndel')
                store_records[0] = record

            vcf_writer.write_record(record)

    finally:
        in_reader.close()
        in_walker.close()
        in_counter.close()
        out_handle.close()
