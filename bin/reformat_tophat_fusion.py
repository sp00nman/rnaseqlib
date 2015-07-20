#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
==========================================
read in potential_fusion.txt & results.txt
==========================================
"""

import sys
import pandas as pd


if __name__ == '__main__':
	
    potential_fusion = sys.argv[1] # potential_fusion.txt
    results_fusion = sys.argv[2] # results.txt
    output_filename = sys.argv[3] # default: merged_data.txt
    

    tophat_potential_fusion = open(potential_fusion, "r")
    
    # except IOError as err:
    #	print('Data file is missing!' + str(err))
    try:
        lines_potential_fusion = tophat_potential_fusion.readlines()
    
    finally:
        tophat_potential_fusion.close()
    
    lines = lines_potential_fusion
    num_lines = range(len(lines))
    for i in num_lines:
	lines[i] = lines[i].strip()
	
        if (i == 0):
	    next
	else:
            if (i % 6)==0:
            #add newline to previous line
                lines[i-1] = lines[i-1] + "\n"
	    else:
	        lines[i] = "\t" + lines[i]
	    i = i + 1

    for i in range(len(lines)):
        
        if ((i+1) % 6)!=0:
            lines[i] = lines[i].replace(' ', '\t')
            
    #print lines

    intermediate_filename = potential_fusion + "_reformat"
    
    reformat_output = open(intermediate_filename, 'wb')
    reformat_output.writelines(lines)
    reformat_output.close()

    # read in potential_fusion_reformat.txt
    df_potential_fusion = pd.read_csv(intermediate_filename, sep="\t",
        names=["sample", "chr", "pos1_breakpoint", "pos2_breakpoint", "orientation",
        "nr_reads_span_fusion", "mate_pairs_support", "mate_pairs_span_fusion",
        "nr_reads_contradict", "nr_reads_left", "nr_reads_right","50bp_contig_pos1_5p",
        "50bp_contig_pos1_3p", "50bp_contig_pos2_5p", "50bp_contig_pos2_3p", "quality_1",
        "quality2", "pos1_gene", "pos1_exon", "pos2_gene", "pos2_exon", "read_depth"])
    
    df_results = pd.read_csv(results_fusion, sep="\t",
        names=["sample", "pos1_gene", "pos1_chr", "pos1_breakpoint", "pos2_gene", "pos2_chr",
        "pos2_breakpoint", "nr_reads_span_fusion", "mate_pairs_support",
        "mate_pairs_span_fusion", "score"])
    
    #if how not specified defaults to inner
    merged_df = pd.merge(df_results, df_potential_fusion, 
            left_on=["sample", "pos1_gene", "pos1_breakpoint", "pos2_gene", "pos2_breakpoint"],
            right_on=["sample", "pos1_gene", "pos1_breakpoint", "pos2_gene", "pos2_breakpoint"],
            sort=True)
    
    if not output_filename:
        output_filename = "merged_data.txt"
    
    merged_df.to_csv(output_filename, sep="\t", index=False)


