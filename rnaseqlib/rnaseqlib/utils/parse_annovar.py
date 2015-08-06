import pandas as pd


def parse_variant(input_file,
                  output_file):
    """
    Filter *.variant_function by
     - only include variants that were flagged as PASS
     - include variants with annotation ncRNA_exonic
     - include variants with annotation splicing
    :return: None
    """
    df = pd.read_csv(
        input_file,
        sep="\t",
        names=["anno", "genesymbol_anno",
               "chrom_anno", "start_anno",
               "end_anno", "ref_anno",
               "alt_anno", "chrom",
               "pos", "id",
               "ref", "alt",
               "qual", "filter",
               "info", "gt",
               "gtnum"])

    df = df[df['filter'] == "PASS"]
    df = df[(df['anno'] == "ncRNA_exonic") | (df['anno'] == "splicing")]
    df.to_csv(output_file, sep="\t", index=0, header=0)

    return None


def parse_exonic(input_file,
                 output_file):
    """
    Only select for variants with flag = PASS
    Introduce a column exonic
    """

    df = pd.read_csv(
        input_file,
        sep="\t",
        names=["line_num", "var",
               "genesymbol_anno","chrom_anno",
               "start_anno", "end_anno",
               "ref_anno", "alt_anno",
               "chrom", "pos",
               "id", "ref",
               "alt", "qual",
               "filter", "info",
               "gt", "gtnum"])

    df = df[df['filter'] == "PASS"]
    df['anno'] = "exonic"
    df = df[['anno', 'var', 'genesymbol_anno','chrom_anno',
             'start_anno', 'end_anno', 'ref_anno', 'alt_anno',
             'chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter',
             'info', 'gt', 'gtnum']]
    df.to_csv(output_file, sep="\t", index=0, header=0)

    return None
