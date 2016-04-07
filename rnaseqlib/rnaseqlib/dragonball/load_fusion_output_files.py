"""
Supports loading fusions to pandas data.frame from tools like
defuse -
tophatfusion
soapfuse
"""

import pandas as pd

DEFUSE_COLUMN = [
    "cluster_id",
    "splitr_sequence",
    "splitr_count", # originally splitr_count
    "splitr_span_pvalue",
    "splitr_pos_pvalue",
    "splitr_min_pvalue",
    "adjacent",
    "altsplice",
    "break_adj_entropy1",
    "break_adj_entropy2",
    "break_adj_entropy_min",
    "breakpoint_homology",
    "breakseqs_estislands_percident",
    "cdna_breakseqs_percident",
    "deletion",
    "est_breakseqs_percident",
    "eversion",
    "exonboundaries",
    "expression1",
    "expression2",
    "gene_A", # originally gene1 (ensembl gene symbol)
    "gene_B", # originally gene2 (ensembl gene symbol)
    "gene_align_strand1",
    "gene_align_strand2",
    "chr_A", # originally gene_chromosome1
    "chr_B", # originally gene_chromosome2
    "gene_end1",
    "gene_end2",
    "gene_location1",
    "gene_location2",
    "gene_name1",
    "gene_name2",
    "gene_start1",
    "gene_start2",
    "gene_strand1",
    "gene_strand2",
    "genome_breakseqs_percident",
    "genomicbreakpoint_A",  # originally genomic_break_pos1
    "genomicbreakpoint_B",  # originally genomic_break_pos2
    "genomic_strand1",
    "genomic_strand2",
    "interchromosomal",
    "interrupted_index1",
    "interrupted_index2",
    "inversion",
    "library_name",
    "max_map_count",
    "max_repeat_proportion",
    "mean_map_count",
    "min_map_count",
    "num_multi_map",
    "num_splice_variants",
    "orf",
    "read_through",
    "repeat_proportion1",
    "repeat_proportion2",
    "span_count",  # originally span_count
    "span_coverage1",
    "span_coverage2",
    "span_coverage_max",
    "span_coverage_min",
    "splice_score",
    "splicing_index1",
    "splicing_index2",
    "score"  # originally probability
]

TOPHATFUSION_COLUMN = [
    "sample",
    "gene_A",  # originally pos1_gene
    "chr_A",  # originally pos1_chr
    "genomicbreakpoint_A",  # originally pos1_breakpoint
    "gene_B",  # originally pos2_gene
    "chr_B",  # originally pos2_chr
    "genomicbreakpoint_B",  # originally pos2_breakpoint
    "splitr_count",  # originally nr_reads_span_fusion_x
    "span_count",  # originally mate_pairs_support_x
    "mate_pairs_span_fusion_x",
    "score",
    "chr",
    "orientation",
    "nr_reads_span_fusion_y",
    "mate_pairs_support_y",
    "mate_pairs_span_fusion_y",
    "nr_reads_contradict",
    "nr_reads_left",
    "nr_reads_right",
    "50bp_contig_pos1_5p",
    "50bp_contig_pos1_3p",
    "50bp_contig_pos2_5p",
    "50bp_contig_pos2_3p",
    "quality_1",
    "quality2",
    "pos1_exon",
    "pos2_exon",
    "read_depth"
]

SOAPFUSE_COLUMN = [
    "gene_A", # originally up_gene
    "up_tran",
    "chr_A", # originally up_chr
    "up_strand",
    "up_Tran_pos",
    "genomicbreakpoint_A", # originally up_Genome_pos
    "up_loc",
    "gene_B", # originally dw_gene
    "dw_tran",
    "chr_B", # originally dw_chr
    "dw_strand",
    "dw_Tran_pos",
    "genomicbreakpoint_B", # originally dw_Genome_pos
    "dw_loc",
    "span_count",  # originally Span_reads_num
    "splitr_count"  # originally Junc_reads_num
]

COLUMN_NAMES = {'defuse': DEFUSE_COLUMN ,
                'tophatfusion': TOPHATFUSION_COLUMN,
                'soapfuse':SOAPFUSE_COLUMN}


def load_fusion(input_file,
                tool):
    """
    Load fusion files from defuse, tophat or soapfuse
    :param input: inputfile
    :param tool: defuse,tophatfusion,sopafuse
    :return: pandas dataframe with fusions
    """

    return pd.read_csv(input_file,
                       sep="\t",
                       header=0,
                       names=COLUMN_NAMES[tool])