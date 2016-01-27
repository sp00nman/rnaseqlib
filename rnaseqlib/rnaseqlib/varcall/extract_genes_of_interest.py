"""
Given a list of genes, extract from table only the genes of interest.
"""

import pandas as pd
from rnaseqlib.utils import tools as ts
from rnaseqlib.utils import convert_gene_ids as cv


VCF_TABLE = \
    ["CHROM",
     "POS",
     "REF",
     "ALT",
     "QUAL",
     "DP",
     "QD",
     "Func.ensGene",
     "Gene_ensGene",
     "ExonicFunc_ensGene",
     "AAChange_ensGene",
     "SIFT_score",
     "Polyphen2_HDIV_score",
     "CADD_phred",
     "GT",
     "AD"]

SELECT_COLUMNS = \
    ["UNIQ_SAMPLE_ID",
     "PATIENT_ID",
     "CHROM",
     "POS",
     "REF",
     "ALT",
     "GENESYMBOL",
     "ENSEMBL_GENEID",
     "ENSEMBL_TRANSCRIPTID",
     "Func.ensGene",
     "ExonicFunc_ensGene",
     "CANONICAL_TRANSCRIPT_EXON_NUM",
     "CANONICAL_TRANSCRIPT_NUCLEOTIDE_CHANGE",
     "CANONICAL_TRANSCRIPT_AS_CHANGE",
     "QUAL",
     "DP",
     "QD",
     "AD",
     "VARIANT_FREQUENCY",
     "SIFT_score",
     "Polyphen2_HDIV_score",
     "CADD_phred"]


def convert_geneids(key_file,
                    ensids):
    """
    For each ensemble gene id convert to genesymbol if available
    :param key_file: dictionary of ensembl gene ids to genesymbol mapping
    :param ensids: list of ensids
    :return:
    """

    genesymbols = []

    for ensid in ensids:
        genesymbol = cv.ensgene2genesymbol(key_file, ensid)
        genesymbols.append(genesymbol)

    return genesymbols


def select_canonical_transcript(
        row,
        canonical_transcripts,
        sep=","):
    """
    Extract from all transcript the canonical transcript
    :param row: for each row in variant_table
    :param canonical_transcripts: a dictionary with key=ensgeneid and
    value=transcriptid
    :param sep: delimiter of that field, "," - is standard
    :return: str with the canonical transcript
    """

    transcripts = row.split(sep)

    # if no canonical transcript is found, take the first...
    cano_trans = transcripts[0]

    for transcript in transcripts:
        if transcript == "UNKNOWN":
            cano_trans = transcript
        else:
            fields = transcript.split(":")

            if fields[0] in canonical_transcripts:
                if fields[1] == canonical_transcripts[fields[0]]:
                 cano_trans = transcript
            else:
                # keep gene_symbol if nothing was found
                cano_trans = transcript
    return cano_trans


def calculate_vf(
        row,
        sep=","):
    """
    Calculate frequency of alternate allele.
    :param row: for each row in variant table
    :param sep: delimiter of REF and ALT read counts; "," - is standard
    :return: float; allele frequency
    """
    return float(row.split(sep)[0])/(float(row.split(sep)[0])+float(row.split(sep)[1]))


def extract_genes(
        project_dir,
        output_file,
        target_list,
        vcf_files,
        conversion_table,
        canonical_transcripts,
        quality,
        depth,
        quality_by_depth,
        exonic_func="nonsynonymous_SNV"):

    """
    Select genes of interest, filter variant and parse output.
    :param project_dir: Project directory
    :param output_file: Name of output file
    :param target_list: List of genes of interest [ENSID\nENSID..]
    :param vcf2table: Dataframe with selected columns.
    :param conversion_table: id conversion table; genesymbol - ensids
    :param exonic_func: nonsynonymous_SNV, synonymous_SNV,stopgain,stoploss,
    non_frameshift_deletion,non_frameshift_insertion,frameshift_deletion,
    frameshift_insertion,unknown
    :param canonical_transcript: id table with mapping of ensids and canonical
    transcript of ensids
    :param quality: QUAL field in VCF file
    :param depth: DP field in VCF file
    :param quality_by_depth: QD field in VCF file
    :return:None
    """

    out_handle = open(output_file, 'a')
    exonic_func_list = exonic_func.split(',')

    for vcf_file in vcf_files:

        #vcf_file[0] has batch attached to it
        uniq_sample_id = vcf_file[0]
        patient_id = vcf_file[1]

        print uniq_sample_id

        # read in vcf2table file into pandas data.frame object
        variant_table = pd.read_csv(
            project_dir + "_"
                + uniq_sample_id
                + ".vcf2table",
            sep="\t",
            names=VCF_TABLE,
            header=True
        )

        #extract only genes of interest
        variant_table = variant_table[variant_table.Gene_ensGene.isin(
            target_list)]

        #extract only exonic functions of interest
        variant_table = variant_table[variant_table.ExonicFunc_ensGene.isin(
            exonic_func_list)]

        #quality filters
        variant_table = variant_table[variant_table['QUAL'] > int(quality)]
        variant_table = variant_table[variant_table['DP'] > int(depth)]
        variant_table = variant_table[variant_table['QD'] > int(quality_by_depth)]

        # add patient information
        variant_table['UNIQ_SAMPLE_ID'] = uniq_sample_id
        variant_table['PATIENT_ID'] = patient_id

        # convert ensids to genesymbols
        genesymbols = convert_geneids(
            cv.read_ensgene_genesymb(conversion_table),
            list(variant_table['Gene_ensGene'])
        )

        # add genesymbol to table
        variant_table['GENESYMBOL'] = genesymbols

        # if dataframe is not empty continue with parsing content
        if not variant_table.empty:

            # calculate allele frequency
            variant_table['VARIANT_FREQUENCY'] = \
                variant_table.apply(lambda row: calculate_vf(
                    row['AD'],
                    sep=","), axis=1
                )

            # only select canonical transcript
            variant_table['CANONICAL_TRANSCRIPT_AAChange'] = \
                variant_table.apply(lambda row: select_canonical_transcript(
                    row['AAChange_ensGene'],
                    canonical_transcripts=canonical_transcripts,
                    sep=","), axis=1
                )

            # split AAChange_ensGene by ':'
            variant_table_split = variant_table[
                'CANONICAL_TRANSCRIPT_AAChange'].apply(lambda row: pd.Series(row.split(':')))
            # rename columns
            variant_table_split.columns = [
                'ENSEMBL_GENEID',
                'ENSEMBL_TRANSCRIPTID',
                'CANONICAL_TRANSCRIPT_EXON_NUM',
                'CANONICAL_TRANSCRIPT_NUCLEOTIDE_CHANGE',
                'CANONICAL_TRANSCRIPT_AS_CHANGE']

            # concatenate the dataframes
            variant_table_concat = pd.concat([variant_table,variant_table_split], axis=1)

            # select specific columns
            variant_table_concat = variant_table_concat[SELECT_COLUMNS]

            variant_table_concat.to_csv(
                out_handle,
                sep="\t",
                index=0,
                header=SELECT_COLUMNS
            )

    out_handle.close()
