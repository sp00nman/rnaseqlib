#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
==============================================================================
Label gene fusions -
==============================================================================
"""

import sys
import optparse


def parse_options(parser):
    parser.add_option("--input_fusion_genes",
                      action="store",
                      type="string",
                      dest="input_fusion_genes_filename",
                      help="""The input file in text tab delimited format
                        containing the fusion genes candidates
                        produced by 'find_fusion_genes.py'. """)

    parser.add_option("--filter_gene_pairs",
                      action="store",
                      type="string",
                      dest="input_filter_gene_pairs_filename",
                      help="""The input text file tab separated format
                           containing gene pairs which are used as
                           filter for labeling (two columns and no header;
                           the order of genes in the gene pairs is ignored).""")

    parser.add_option("--filter_genes",
                      action="store",
                      type="string",
                      dest="input_filter_genes_filename",
                      help="""The input text file format containing genes which
                           are used as filter for labeling.""")

    parser.add_option("--label",
                      action="store",
                      type="string",
                      dest="label",
                      help="""Label used to mark the candidate fusion genes
                            which are founf in the filter.""")

    parser.add_option("--output_fusion_genes",
                      action="store",
                      type="string",
                      dest="output_fusion_genes_filename",
                      help="""The output text tab-separated file containing
                           the candidate fusion genes which are found in the
                           filter. The format is as the input file and sorted
                           by counts column.""")

    parser.add_option("--fusion_detection_algorithm",
                      action="store",
                      type="string",
                      dest="fusion_detection_algorithm",
                      help="""What is the input fusion det. tool?.""")

    parser.add_option("--min_dist_gene_gene",
                      action="store",
                      type="int",
                      dest="input_min_dist_gene_gene",
                      help="""Labels pairs of genes below a given threshold.""")

    parser.add_option("--min_dist_gene_gene_database",
                      action="store",
                      type="string",
                      dest="input_min_dist_gene_gene_database_filename",
                      help="""Database with exons position on chromosomes which is used
                           to extract the gene positions, e.g. 'more_exons_ensembl.txt'.
                           This is needed only when '--min_distance_gene_gene' is used.""")

    (options, args) = parser.parse_args()

    print "Input parameters"
    print "[1] input_fusion_genes_filename: ", options.input_fusion_genes_filename
    print "[2] input_filter_gene_pairs_filename: ", options.input_filter_gene_pairs_filename
    print "[3] output_fusion_genes_filename:", options.output_fusion_genes_filename
    print "[4] label:", options.label
    print "[5] fusion_detection_algorithm: ", options.fusion_detection_algorithm, "\n"

    return options


def validate_options(options, parser):
    """
    Validate user input parameters
    """
    if not ((options.input_fusion_genes_filename and options.output_fusion_genes_filename) and options.label):
        parser.print_help()
        parser.error("One of the options has not been specified.")
        sys.exit(1)


def read_gene_pairs(options):
    """
    gene pairs: eg. ENSG....\tENSG....
    """
    gene_pairs = set()
    if options.input_filter_gene_pairs_filename:
        print "Reading,...", options.input_filter_gene_pairs_filename
        gene_pairs = set(['\t'.join(sorted(line.rstrip('\r\n').split('\t')[:2]))
                          for line in file(options.input_filter_gene_pairs_filename, 'r') if line.rstrip('\r\n')])
    return gene_pairs


def read_no_proteins(options):
    """
    TODO: improve description of function
    filter genes for non-coding ?
    """
    no_proteins = set()
    if options.input_filter_genes_filename:
        print "Reading...", options.input_filter_genes_filename
        no_proteins = set([line.rstrip('\r\n')
                           for line in file(options.input_filter_genes_filename, 'r') if line.rstrip('\r\n')])
    return no_proteins


def read_exon_database(options):
    """
    ensembl exon database
    ensembl_peptide_id             0
    ensembl_gene_id                1
    ensembl_transcript_id          2
    ensembl_exon_id                3
    exon_chrom_start               4
    exon_chrom_end                 5
    rank                           6
    start_position                 7
    end_position                   8
    transcript_start               9
    transcript_end                 10
    strand                         11
    chromosome_name                12
    cds_start                      13
    cds_end                        14
    5_utr_start                    15
    5_utr_end                      16
    3_utr_start                    17
    3_utr_end                      18
    """
    print "Processing the exons database...", options.input_min_dist_gene_gene_database_filename
    genes = {}
    if options.input_min_dist_gene_gene_database_filename:
        exons = [line.rstrip('\r\n').split('\t')
                 for line in file(options.input_min_dist_gene_gene_database_filename, 'r').readlines()
                 if line.rstrip('\r\n')]

        exons = [(line[1],
                  line[7],
                  line[8],
                  line[11],
                  line[12])
                 for line in exons]

        for line in exons:
            gn = line[0]
            gs = int(line[1])
            ge = int(line[2])
            st = int(line[3])
            ch = line[4]
            if gs > ge:
                (gs, ge) = (ge, gs)
            if gn not in genes:
                genes[gn] = {'start': gs,
                             'end': ge,
                             'strand': st,
                             'chrom': ch}

    return genes


def read_fusion_input(options):
    """
    Read in fusion detection file input 21 22: gene1 bzw gene ?
    """
    print "Reading...", options.input_fusion_genes_filename
    data = [line.rstrip('\r\n').split('\t')
            for line in file(options.input_fusion_genes_filename, 'r').readlines() if line.rstrip('\r\n')]

    return data


def read_ensembl_id(path2files):
    """
    old file : hg19_ensgene_simple_4labelfusions_sortuniq_woENST.txt
    new file : cut_merged_gene2ensembl_geneinfo_human
    Entrez_GeneID gene_symbol	ensemblID
    1	A1BG	ENSG00000121410
    """
    ensembl_id_genesymbol_f = path2files + "cut_merged_gene2ensembl_geneinfo_human"
    print "Reading...", ensembl_id_genesymbol_f
    ensembl_id_genesymbol = {}
    ensembl_id_genesymbol_data = [line.rstrip('\r\n').split('\t')
                                  for line in file(ensembl_id_genesymbol_f, 'r').readlines() if line.rstrip('\r\n')]

    for line in ensembl_id_genesymbol_data:
        ensembl_id = line[2]
        genesymbol = line[1]
        ensembl_id_genesymbol[genesymbol] = ensembl_id

    return ensembl_id_genesymbol


def read_missing_data(path2files):
    """
    fill dictionary
    key = ensemblID | value = genesymobls
    careful file contains 57389 entries
    genesymbols are not unique, after read in 50480 entries remain
    """
    missing_data_f = path2files + "mod_diff.txt"
    print "Reading...", missing_data_f
    missing_data = {}
    missing_data_input = [line.rstrip('\r\n').split('\t')
                          for line in file(missing_data_f, 'r').readlines() if line.rstrip('\r\n')]

    for line in missing_data_input:
        ensembl_id_miss = line[0]
        genesymbol_miss = line[1]
        missing_data[genesymbol_miss] = ensembl_id_miss

    return missing_data


def read_refseqid(path2files):
    """
    Read in refseqid_ensemblid_f
    filename = homo_sapiens_core_73_37_NM_NR_numbers_ensembl_stabl_id.txt
    key = ensemblid | value = refseqid
    """
    refseqid_ensemblid_f = path2files + "homo_sapiens_core_73_37_NM_NR_numbers_ensembl_stabl_id.txt"
    print "Reading...", refseqid_ensemblid_f
    refseqid_ensemblid = {}
    refseqid_ensemblid_input = [line.rstrip('\r\n').split('\t')
                                for line in file(refseqid_ensemblid_f, 'r').readlines() if line.rstrip('\r\n')]

    for line in refseqid_ensemblid_input:
        ensemblid = line[2]
        refseqid = line[1]
        refseqid_ensemblid[refseqid] = ensemblid

    return refseqid_ensemblid


def gene_symbol2ensembl_id(ensemblid_genesymbol, missing_data, gene_symbol):
    """
    fuction for converting gene_symbols to ensembl_IDs
    careful: LOC541471 changed
    """
    if gene_symbol in ensemblid_genesymbol:
        ensemblid = ensemblid_genesymbol[gene_symbol]
    elif gene_symbol in missing_data:
        ensemblid = missing_data[gene_symbol]
    else:
        ensemblid = gene_symbol  # keep gene_symbol if nothing was found

    return ensemblid


def refseqid2ensemblid(refseqid_ensemblid, refseqid):
    if refseqid in refseqid_ensemblid:
        ensemblid = refseqid_ensemblid[refseqid]
    else:
        ensemblid = refseqid  # keep gene_symbol if nothing was found
    return ensemblid


def get_map_genes(line, options, col_num_gene1, col_num_gene2):
    if options.fusion_detection_algorithm == 'defuse':
        a = line[col_num_gene1[options.fusion_detection_algorithm]]
        b = line[col_num_gene2[options.fusion_detection_algorithm]]

    elif options.fusion_detection_algorithm == 'fusionmap':
        split_gene1 = line[col_num_gene1[options.fusion_detection_algorithm]].split(",")
        split_gene2 = line[col_num_gene2[options.fusion_detection_algorithm]].split(",")

        #TODO: make a for loop here
        if type(split_gene1) is list:
            split_refid = split_gene1[0].split("_")[0] + "_" + split_gene1[0].split("_")[1]
            a = refseqid2ensemblid(split_refid)

        else:
            a = refseqid2ensemblid(line[col_num_gene1[options.fusion_detection_algorithm]])

        if type(split_gene2) is list:
            split_refid = split_gene2[0].split("_")[0] + "_" + split_gene2[0].split("_")[1]
            b = refseqid2ensemblid(split_refid)

        else:
            b = refseqid2ensemblid(line[col_num_gene1[options.fusion_detection_algorithm]])

    else:
        a = gene_symbol2ensembl_id(line[col_num_gene1[options.fusion_detection_algorithm]])
        b = gene_symbol2ensembl_id(line[col_num_gene2[options.fusion_detection_algorithm]])
    return a, b


def min_dis_gene(data, genes, options, label_col,
                 col_num_gene1, col_num_gene2):
    """
    deal with the distance between genes and label accordingly
    """
    temp = []
    for line in data:

        (a, b) = get_map_genes(line, options, col_num_gene1, col_num_gene2)

        if (a in genes and
                b in genes and
                    genes[a]["chrom"] == genes[b]["chrom"] and
                    genes[a]["strand"] == genes[b]["strand"] and
                    min([abs(genes[a]["start"] - genes[b]["start"]),
                         abs(genes[a]["start"] - genes[b]["end"]),
                         abs(genes[a]["end"] - genes[b]["start"]),
                         abs(genes[a]["end"] - genes[b]["end"])]) <= options.input_min_dist_gene_gene):

            if label_col:
                temp.append(line + [options.label])
            else:
                if line[-1]:
                    temp.append(line[:-1] + [','.join([line[-1], options.label])])
                else:
                    temp.append(line[:-1] + [options.label])
        else:
            if label_col:
                temp.append(line + [''])
            else:
                temp.append(line)
    return temp


def filter_gene_pairs(data, gene_pairs, no_proteins,
                      col_num_gene1, col_num_gene2,
                      label_col, options):
    """
    for gene pairs
    """
    temp = []
    for line in data:
        (a, b) = get_map_genes(line, options, col_num_gene1, col_num_gene2)
        g = '\t'.join(sorted([a, b]))
        print g

        if (g in gene_pairs) or (a in no_proteins) or (b in no_proteins):
            print g
            if label_col:
                temp.append(line + [options.label])
            else:
                if line[-1]:
                    temp.append(line[:-1] + [','.join([line[-1], options.label])])
                else:
                    temp.append(line[:-1] + [options.label])
        else:
            if label_col:
                temp.append(line + [''])
            else:
                temp.append(line)
    return temp


def main(argv=None):
    path2files = "/home/illumina/Fiorella/annotation/"
    col_num_total = {'defuse': 65, 'fusionmap': 21, 'tophat-fusion': 27}
    col_num_gene1 = {'defuse': 20, 'fusionmap': 10, 'tophat-fusion': 1}
    col_num_gene2 = {'defuse': 21, 'fusionmap': 14, 'tophat-fusion': 3}

    usage = "%prog [options]"
    description = """Short: Labels fusions"""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,
                                   description=description,
                                   version=version)

    options = parse_options(parser)
    validate_options(options, parser)

    data = read_fusion_input(options)
    label_col = False
    header = data.pop(0)
    if len(header) == col_num_total[options.fusion_detection_algorithm]:
        label_col = True
        header.append('Fusion_description')

    read_ensembl_id(path2files)
    missing_data = read_missing_data(path2files)
    read_refseqid(path2files)
    no_proteins = read_no_proteins(options)

    if options.input_min_dist_gene_gene:
        temp = min_dis_gene(data, label_col, options,
                            col_num_gene1, col_num_gene2)
        genes = read_exon_database(options)
    else:
        gene_pairs = read_gene_pairs(options)
        temp = filter_gene_pairs(data, gene_pairs, no_proteins,
                                 label_col, options, col_num_gene1, col_num_gene2)

    data = sorted(temp)
    data.insert(0, header)
    file(options.output_fusion_genes_filename, 'w').writelines(['\t'.join(line) + '\n'
                                                                for line in data])


if __name__ == "__main__":
    sys.exit(main())
