"""
Filter vcf (vcf file needs to be ANNOVAR annotated with the following
databases: ensGene,snp138NonFlagged,1000g2014oct_all
"""

import vcf


def filter_vcf(input_file,
               file_prefix,
               file_suffix):
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

    Statistics: Create table with the following information:
    - COUNT_PASS
    - COUNT_FS
    - COUNT_QD
    - COUNT_SNPCLUSTER
    - COUNT_HRUN
    - COUNT NHRUN
    - COUNT_TOTAL
    - COUNT_INDEL
    - COUNT_SNV
    - COUNT_EXONIC
    - COUNT_INTRONIC
    - COUNT_SPLICING
    - COUNT_UTR
    - COUNT_EXONIC_SPLICING
    - COUNT_NON_CODING
    - COUNT_RARE_VARIANTS
    - COUNT_SYNONYMOUS

    :param input_file: Input VCF file
    :param file_suffix: suffix VCF file
    :return:None
    """

    vcf_file = open(input_file, 'r')
    stats = open(file_prefix + ".stats_variants.txt", 'w')
    # separate tables
    all = open(file_prefix + ".all." + file_suffix, 'w')
    exonic = open(file_prefix + ".exonic." + file_suffix, 'w')
    splicing = open(file_prefix + ".splicing." + file_suffix, 'w')
    UTR = open(file_prefix + ".UTR." + file_suffix, 'w')
    exonic_splicing = open(file_prefix + ".exonic_splicing." + file_suffix, 'w')
    ncRNA_exonic = open(file_prefix + ".ncRNA_exonic." + file_suffix, 'w')
    ncRNA_splicing = open(file_prefix + ".ncRNA_splicing." + file_suffix, 'w')

    try:
        vcf_reader = vcf.Reader(vcf_file)
        vcf_all = vcf.Writer(all, vcf_reader)
        vcf_exonic = vcf.Writer(exonic, vcf_reader)
        vcf_splicing = vcf.Writer(splicing, vcf_reader)
        vcf_UTR = vcf.Writer(UTR, vcf_reader)
        vcf_exonic_splicing = vcf.Writer(exonic_splicing, vcf_reader)
        vcf_ncRNA_exonic = vcf.Writer(ncRNA_exonic, vcf_reader)
        vcf_ncRNA_splicing = vcf.Writer(ncRNA_splicing, vcf_reader)


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

            # filter variants
            # for some odd reason is [None] doesn't work, is that
            # an empty list ?
            if not record.FILTER \
                    and record.INFO['snp142Common'] == [None] \
                    and record.INFO['1000g2015feb_all'] < 0.01 \
                    and record.INFO['esp5400_all'] < 0.01 \
                    and record.INFO['esp6500siv2_all'] < 0.01 \
                    and record.INFO['ExonicFunc.ensGene'][0] \
                            != "synonymous_SNV":

                COUNT_RARE_VARIANTS += 1
                vcf_all.write_record(record)

                for info in record.INFO['Func.ensGene']:
                    if info == "exonic":
                        vcf_exonic.write_record(record)
                    if info == "splicing":
                        vcf_splicing.write_record(record)
                    if info == "UTR3" or info == "UTR5":
                        vcf_UTR.write_record(record)
                    if info == "exonic;splicing":
                        vcf_exonic_splicing.write_record(record)
                    if info == "ncRNA_exonic" :
                        vcf_ncRNA_exonic.write_record(record)
                    if info == "ncRNA_splicing":
                        vcf_ncRNA_splicing.write_record(record)

        # write statistics to file
        text = "SAMPLE_NAME"                    + "\t" \
                + "COUNT_PASS"                  + "\t" \
                + "COUNT_FS"                    + "\t" \
                + "COUNT_QD"                    + "\t" \
                + "COUNT_SNPCLUSTER"            + "\t" \
                + "COUNT_LOWQUAL"               + "\t" \
                + "COUNT_HRUN"                  + "\t" \
                + "COUNT_NHRUN"                 + "\t" \
                + "COUNT_NINDEL"                + "\t" \
                + "COUNT_TOTAL"                 + "\t" \
                + "COUNT_INDEL"                 + "\t" \
                + "COUNT_SNV"                   + "\t" \
                + "COUNT_EXONIC"                + "\t" \
                + "COUNT_INTRONIC"              + "\t" \
                + "COUNT_SPLICING"              + "\t" \
                + "COUNT_UTR"                   + "\t" \
                + "COUNT_EXONIC_SPLICING"       + "\t" \
                + "COUNT_NON_CODING"            + "\t" \
                + "COUNT_RARE_VARIANTS"         + "\t" \
                + "COUNT_SYNONYMOUS"            + "\n" \
                + file_prefix                   + "\t" \
                + str(COUNT_PASS)               + "\t" \
                + str(COUNT_FS)                 + "\t" \
                + str(COUNT_QD)                 + "\t" \
                + str(COUNT_SNPCLUSTER)         + "\t" \
                + str(COUNT_LOWQUAL)            + "\t" \
                + str(COUNT_HRUN)               + "\t" \
                + str(COUNT_NHRUN)              + "\t" \
                + str(COUNT_NINDEL)             + "\t" \
                + str(COUNT_TOTAL)              + "\t" \
                + str(COUNT_INDEL)              + "\t" \
                + str(COUNT_SNV)                + "\t" \
                + str(COUNT_EXONIC)             + "\t" \
                + str(COUNT_INTRONIC)           + "\t" \
                + str(COUNT_SPLICING)           + "\t" \
                + str(COUNT_UTR)                + "\t" \
                + str(COUNT_EXONIC_SPLICING)    + "\t" \
                + str(COUNT_NON_CODING)         + "\t" \
                + str(COUNT_RARE_VARIANTS)      + "\t" \
                + str(COUNT_SYNONYMOUS)

        stats.write(text)

    finally:
        vcf_file.close()
        all.close()
        exonic.close()
        splicing.close()
        exonic_splicing.close()
        UTR.close()
        stats.close()

    return None