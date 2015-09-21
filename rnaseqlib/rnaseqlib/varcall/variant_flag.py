"""
Flag variants that are listed in file with the following format:
## generic variant information
CHROM eg. 1
POS eg. 44134887
REF eg. T
ALT eg. C
Gene_ensGene eg. ENSG00000066135
GENESYMBOL eg. KDM4A

## followed by cohort data, all ";" separated
UNIQ_SAMPLE_ID_COHORT eg. CEMM_D01_CC2;CEMM_05_co_CC1;
QUAL_COHORT eg. 606.77;246.77;
DP_COHORT eg. 58;27;
QD_COHORT eg. 10.46;9.14;
ExonicFunc_ensGene_COHORT eg. nonsynonymous_SNV;nonsynonymous_SNV;
GT_COHORT T/C;T/C;
AD_COHORT 32,26;15,12;
"""

import vcf
import pandas as pd


def annotate_variants(anno_file,
                      input_file,
                      output_file,
                      min_num_occurrence,
                      min_vaf):

    in_handle = open(input_file, 'r')
    out_handle = open(output_file, 'w')

    d = {}
    # load anno_file
    variant_table = pd.read_csv(
        anno_file,
        sep="\t",
        names=["CHROM",
               "POS",
               "REF",
               "ALT",
               "Gene_ensGene",
               "GENESYMBOL",
               "UNIQ_SAMPLE_ID_COHORT",
               "QUAL_COHORT",
               "DP_COHORT",
               "QD_COHORT",
               "ExonicFunc_ensGene_COHORT",
               "GT_COHORT",
               "AD_COHORT"],
        header=True
    )

    # read in vcf to nested hash table
    for index, row in variant_table.iterrows():

        # filter variants
        vaflist = []
        num_occurrences = len(str(variant_table.loc[index]['UNIQ_SAMPLE_ID_COHORT']).split(";"))-1

        for AD in variant_table.loc[index]['AD_COHORT'].split(";"):
            if AD != '':
                allele_count = AD.split(",")
                REF = int(allele_count[0])
                ALT = int(allele_count[1])
                vaflist.append(ALT/(ALT+REF))
        max_vaflist = max(vaflist)

        if num_occurrences >= min_num_occurrence and max_vaflist >= min_vaf:
            # fill the hash
            key_variant = str(row['CHROM']) + ":" + str(row['POS']) + ";" \
                         + row['REF'] + ";" + row['ALT']
            value_variant = {'QUAL_COHORT':row['QUAL_COHORT'],
                             'DP_COHORT':row['DP_COHORT'],
                             'QD_COHORT':row['QD_COHORT'],
                             'Gene_ensGene':row['Gene_ensGene'],
                             'GENESYMBOL':row['GENESYMBOL'],
                             'UNIQ_SAMPLE_ID_COHORT':row['UNIQ_SAMPLE_ID_COHORT'],
                             'ExonicFunc_ensGene_COHORT':row['ExonicFunc_ensGene_COHORT'],
                             'GT_COHORT':row['GT_COHORT'],
                             'AD_COHORT':row['AD_COHORT']}

            d[key_variant] = value_variant

    try:
        vcf_reader = vcf.Reader(in_handle)
        vcf_writer = vcf.Writer(out_handle, vcf_reader)

        for record in vcf_reader:
            # record.ALT returns a list,...weird hack needed to
            # get the correct type for comparison
            REF_allele = str(record.alleles[0])
            ALT_allele = str(record.alleles[1])
            key = str(record.CHROM) + ":" + str(record.POS) \
                  + ";" + str(REF_allele) + ";" + str(ALT_allele)
            if key in d:
                # change CC to other abbreviation, easier to grep...
                record.FILTER.append('CC_ALL')

            vcf_writer.write_record(record)

    finally:
        in_handle.close()
        out_handle.close()