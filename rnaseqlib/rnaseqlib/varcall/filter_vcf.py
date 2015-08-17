"""
Filter vcf (vcf file needs to be ANNOVAR annotated with the following
databases: ensGene,snp138NonFlagged,1000g2014oct_all
"""

import vcf


def filter_vcf(input_file,
               output_file):
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
    - COUNT TOTAL
    - COUNT_INDEL
    - COUNT_SNV
    - COUNT_EXONIC
    - COUNT_INTRONIC
    - COUNT_SPLICING
    - COUNT_UTR
    - COUNT EXONIC_SPLICING
    - COUNT NON_CODING
    - COUNT RARE_VARIANTS
    - COUNT_SYNONYMOUS

    :param input_file: Input VCF file
    :param output_file: suffix VCF file
    :return:None
    """

    vcf_file = open(input_file, 'r')
    out_file = open(output_file, 'w')

    try:
        vcf_reader = vcf.Reader(vcf_file)
        vcf_writer = vcf.Writer(out_file, vcf_reader)

        # define variables:
        COUNT_PASS = 0
        COUNT_FS = 0
        COUNT_QD = 0
        COUNT_SNPCLUSTER = 0
        COUNT_LOWQUAL = 0
        COUNT_HRUN = 0
        COUNT_NHRUN = 0
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

            if record.is_indel:
                COUNT_INDEL += 1

            if not record.is_indel:
                COUNT_SNV += 1

            if record.INFO['Func.ensGene'] is "exonic":
                COUNT_EXONIC += 1

            if record.INFO['Func.ensGene'] is "splicing":
                COUNT_SPLICING += 1

            if record.INFO['Func.ensGene'] is "UTR":
                COUNT_UTR += 1

            if record.INFO['ExonicFunc.ensGene'] is "synonymous_SNV":
                COUNT_SYNONYMOUS += 1

            # filter variants
            if not record.FILTER:
                if not record.INFO['snp138NonFlagged']:
                    if record.INFO['1000g2014oct_all'] < 0.1:
                        if record.INFO['ExonicFunc.ensGene'] \
                                is not "synonymous_SNV":
                            COUNT_RARE_VARIANTS += 1
                            vcf_writer.write_record(record)

        print COUNT_PASS, COUNT_FS, COUNT_QD, COUNT_SNPCLUSTER, COUNT_LOWQUAL,\
            COUNT_HRUN, COUNT_NHRUN, COUNT_TOTAL, COUNT_INDEL, COUNT_SNV, \
            COUNT_EXONIC, COUNT_INTRONIC, COUNT_SPLICING, COUNT_UTR, \
            COUNT_EXONIC_SPLICING, COUNT_NON_CODING, COUNT_RARE_VARIANTS, \
            COUNT_SYNONYMOUS

    finally:
        vcf_file.close()
        out_file.close()


    return None