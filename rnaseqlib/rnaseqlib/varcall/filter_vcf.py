"""
Filter VCF file:
VCF file needs to be ANNOVAR annotated with the following
databases: ensGene,snp138NonFlagged,1000g2014oct_all
"""

import vcf


def filter_vcf(input_file,
               file_prefix,
               file_suffix,
               mode):
    """
    Filter vcf with the following criteria:
    - FILTER: only keep 'PASS' annotation
    - INFO: snp138Flagged has rs-number
    - INFO: 1000g2014oct_all has AF>0.1
    - INFO: ExonicFunc.ensGene = 'synonymous SNV'

    Create [separate] tables for
    - exonic
    - splicing
    - UTR
    - exonic;splicing
    - ncRNA_exonic
    - ncRNA_splicing

    :param input_file: Input VCF file
    :param file_prefix: prefix of VCF file
    :param file_suffix: suffix VCF file
    :param mode: SOMATIC(filter for somatic variants); CLONXCHR(filter for
    polymorphisms in the X chromosome for assessment fo clonality; ALL(filter
    for all)
    :return:None
    """

    # define variables for statistic
    # would look much nice if programmed object oriented...
    COUNT_PASS = 0
    COUNT_FS = 0
    COUNT_QD = 0
    COUNT_SNPCLUSTER = 0
    COUNT_LOWQUAL = 0
    COUNT_HRUN = 0
    COUNT_NHRUN = 0
    COUNT_NINDEL = 0
    COUNT_TOTAL= 0
    COUNT_INDEL = 0
    COUNT_SNV = 0
    COUNT_EXONIC = 0
    COUNT_INTRONIC = 0
    COUNT_SPLICING = 0
    COUNT_UTR = 0
    COUNT_EXONIC_SPLICING = 0
    COUNT_NON_CODING = 0
    COUNT_RARE_VARIANTS = 0
    COUNT_SYNONYMOUS = 0

    # input file
    vcf_file = open(
        input_file, 'r')

    # output files
    stats = open(
        file_prefix
        + ".stats_variants.txt", 'w')

    if mode == "SOMATIC" or mode == "ALL":
        UTR = open(
            file_prefix
            + ".UTR."
            + file_suffix, 'w')
        exonic_splicing = open(
            file_prefix
            + ".exonic_splicing."
            + file_suffix, 'w')
        ncRNA_exonic = open(
            file_prefix
            + ".ncRNA_exonic."
            + file_suffix, 'w')
        ncRNA_splicing = open(
            file_prefix
            + ".ncRNA_splicing."
            + file_suffix, 'w')
        all = open(
            file_prefix
            + ".all."
            + file_suffix, 'w')

    elif mode == "CLONXCHR" or mode == "ALL":

        clonxchr = open(
            file_prefix
            + ".clonxchr."
            + file_suffix, 'w')

    else:
        print "Unknown mode selected. Options are SOMATIC, ALL, CLONXCHR."

    try:
        vcf_reader = vcf.Reader(vcf_file)

        if mode == "SOMATIC" or mode == "ALL":
            vcf_UTR = vcf.Writer(UTR, vcf_reader)
            vcf_exonic_splicing = vcf.Writer(exonic_splicing, vcf_reader)
            vcf_ncRNA_exonic = vcf.Writer(ncRNA_exonic, vcf_reader)
            vcf_ncRNA_splicing = vcf.Writer(ncRNA_splicing, vcf_reader)
            vcf_all = vcf.Writer(all, vcf_reader)

        elif mode == "CLONXCHR" or mode == "ALL":
            vcf_clonxchr = vcf.Writer(clonxchr, vcf_reader)

        else:
            print "Unknown mode selected. Options are SOMATIC, ALL, CLONXCHR."

        for record in vcf_reader:

            COUNT_TOTAL += 1

            if not record.FILTER:
                COUNT_PASS += 1

            for filter in record.FILTER:
                if filter == 'FS':
                    COUNT_FS += 1
                if filter == 'QD':
                    COUNT_QD += 1
                if filter == 'SnpCluster':
                    COUNT_SNPCLUSTER += 1
                if filter == 'LowQual':
                    COUNT_LOWQUAL += 1
                if filter == 'HRun':
                    COUNT_HRUN += 1
                if filter == 'nHRun':
                    COUNT_NHRUN += 1
                if filter == 'nIndel':
                    COUNT_NINDEL += 1

            if record.is_indel:
                COUNT_INDEL += 1

            if not record.is_indel:
                COUNT_SNV += 1

            for info in record.INFO['Func.ensGene']:
                if info == "intronic":
                    COUNT_INTRONIC += 1
                if info == "exonic":
                    COUNT_EXONIC += 1
                if info == "splicing":
                    COUNT_SPLICING += 1
                if info == "UTR3" or info == "UTR5":
                    COUNT_UTR += 1
                if info == "exonic;splicing":
                    COUNT_EXONIC_SPLICING += 1
                if info == "ncRNA_exonic" or info == "ncRNA_splicing":
                    COUNT_NON_CODING += 1

            for info in record.INFO['ExonicFunc.ensGene']:
                if info == "synonymous_SNV":
                    COUNT_SYNONYMOUS += 1

            if mode == "SOMATIC" or mode == "ALL":

                # for some odd reason is [None] doesn't work, is that
                # an empty list ?
                # if not record.FILTER keeps only PASS variants
                if not record.FILTER \
                        and record.INFO['snp142Common'] == [None] \
                        and record.INFO['1000g2015feb_all'] < 0.01 \
                        and record.INFO['esp5400_all'] < 0.01 \
                        and record.INFO['esp6500siv2_all'] < 0.01 \
                        and record.INFO['ExonicFunc.ensGene'][0] \
                                != "synonymous_SNV":

                    # split records into separate files
                    for info in record.INFO['Func.ensGene']:
                        if info == "exonic" or info == "splicing" \
                                or info == "exonic;splicing":
                            vcf_exonic_splicing.write_record(record)
                        if info == "UTR3" or info == "UTR5":
                            vcf_UTR.write_record(record)
                        if info == "ncRNA_exonic" :
                            vcf_ncRNA_exonic.write_record(record)
                        if info == "ncRNA_splicing":
                            vcf_ncRNA_splicing.write_record(record)

                    # write all records for SOMATIC
                    vcf_all.write_record(record)

                    COUNT_RARE_VARIANTS += 1

            elif mode == "CLONXCHR" or mode == "ALL":
                if not record.FILTER and record.is_snp and record.CHROM == "X":
                    # sorting exonic regions will be done
                    for info in record.INFO['Func.ensGene']:
                        if info == "exonic" \
                                or info == "splicing" \
                                or info =="exonic;splicing" \
                                or info == "UTR3" \
                                or info=="UTR5" \
                                or info == "ncRNA_exonic" \
                                or info == "ncRNA_splicing":
                            vcf_clonxchr.write_record(record)

            else:
                print "Unknown mode selected. Options are SOMATIC, ALL, CLONXCHR."

        # write statistics to file
        STATS_COLUMNNAMES = [
            "SAMPLE_NAME",
            "COUNT_PASS",
            "COUNT_FS",
            "COUNT_QD" ,
            "COUNT_SNPCLUSTER",
            "COUNT_LOWQUAL" ,
            "COUNT_HRUN",
            "COUNT_NHRUN",
            "COUNT_NINDEL",
            "COUNT_SNV",
            "COUNT_EXONIC",
            "COUNT_INTRONIC",
            "COUNT_SPLICING",
            "COUNT_UTR",
            "COUNT_EXONIC_SPLICING",
            "COUNT_NON_CODING",
            "COUNT_RARE_VARIANTS",
            "COUNT_SYNONYMOUS"
        ]

        STATS_VARIABLES = [
            file_prefix,
            str(COUNT_PASS),
            str(COUNT_FS),
            str(COUNT_QD),
            str(COUNT_SNPCLUSTER),
            str(COUNT_LOWQUAL),
            str(COUNT_HRUN),
            str(COUNT_NHRUN),
            str(COUNT_NINDEL),
            str(COUNT_TOTAL),
            str(COUNT_INDEL),
            str(COUNT_SNV),
            str(COUNT_EXONIC),
            str(COUNT_INTRONIC),
            str(COUNT_SPLICING),
            str(COUNT_UTR),
            str(COUNT_EXONIC_SPLICING),
            str(COUNT_NON_CODING),
            str(COUNT_RARE_VARIANTS),
            str(COUNT_SYNONYMOUS)
        ]

        STATS = '\t'.join(STATS_COLUMNNAMES) \
                + "\n" \
                + '\t'.join(STATS_VARIABLES)

        stats.write(STATS)

    finally:
        vcf_file.close()
        stats.close()

        if mode == "SOMATIC" or mode == "ALL":
            all.close()
            exonic_splicing.close()
            UTR.close()
            vcf_ncRNA_exonic.close()
            vcf_ncRNA_splicing.close()

        elif mode == "CLONXCHR" or mode == "ALL":
            vcf_clonxchr.close()

        else:
            print "Unknown mode selected. Options are SOMATIC, ALL, CLONXCHR."

    return None