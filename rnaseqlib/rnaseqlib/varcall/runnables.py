"""
Collection of functions that execute external software or UNIX commands.
"""


def rnaseq_align(star_genome,
                 read1,
                 read2,
                 num_cpus,
                 outfile_prefix):
    """
    source: http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq
    https://github.com/alexdobin/STAR
    Read alignment: STAR 2-pass method which was described in a recent
    publication (see page 43 of the Supplemental text of the Paer G Engstroem
    et al. paper referenced below for full protocol details -- we used the
    suggested protocol with the default parameters). In brief, in the
    STAR 2-pass approach, splice junctions detected in a first alignment
    run are used to guide the final alignment.
    """

    cmd_align = "STAR " \
                "--genomeDir %s " \
                "--readFilesIn %s %s " \
                "--runThreadN %s " \
                "--genomeLoad NoSharedMemory " \
                "--outFilterIntronMotifs RemoveNoncanonical " \
                "outSAMtype BAM" \
                "--outFileNamePrefix %s" % (star_genome, read1,
                                            read2, num_cpus,
                                            outfile_prefix)
    return cmd_align


def star_index(novel_ref,
               genome,
               firstroundalignment,
               sjdb_overhang,
               num_cpus,
               outfile_prefix):
    """
    For the 2-pass STAR, a new index is then created using splice junction
    information contained in the file SJ.out.tab from the first pass.
    """
    cmd_star_index = "STAR " \
                     "--runMode genomeGenerate " \
                     "--genomeDir %s " \
                     "--genomeFastaFiles %s.fa " \
                     "--sjdbFileChrStartEnd %s " \
                     "--sjdbOverhang %s " \
                     "--runThreadN %s " \
                     "--outFileNamePrefix %s" % (novel_ref,
                                                 genome,
                                                 firstroundalignment,
                                                 sjdb_overhang,
                                                 num_cpus,
                                                 outfile_prefix)
    return cmd_star_index


def extract(input_file,
            region,
            output_file):
    """
    Extract a region from a BAM file.
    """

    cmd_extract = "samtools view -b -h %s %s >%s " % (input_file,
                                                      region,
                                                      output_file)
    return cmd_extract


def reorder_sam(inbamfile,
                outbamfile,
                genome_path):
    """
    Reorder BAM file.
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file
    :param genomes: path to genome
    :return: Command to be executed; type str
    """

    cmd_reorder = "java -Xmx6g -jar $NGS_PICARD/ReorderSam.jar " \
                  "INPUT=%s " \
                  "OUTPUT=%s " \
                  "REFERENCE=%s.fa" % (inbamfile,
                                    outbamfile,
                                    genome_path)
    return cmd_reorder


def sort_bam(inbamfile,
             outbamfile,
             sort_order="coordinate"):
    """
    Sort BAM file by variable
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file (BAM formatted and sorted)
    :param sort_order: sort by (default: coordinate)
    :return: Command to be executed; type str
    """

    cmd_sort = "java -Xmx6g -jar $NGS_PICARD/SortSam.jar " \
               "INPUT=%s " \
               "OUTPUT=%s " \
               "SORT_ORDER=%s" % (inbamfile,
                                  outbamfile,
                                  sort_order)
    return cmd_sort


def replace_readgroups(input_file,
                       output_file,
                       project_name):
    """
    Replace readgroups (with some dummy variables.)
    """

    cmd_replace = "java -jar $NGS_PICARD/AddOrReplaceReadGroups.jar " \
                  "I=%s " \
                  "O=%s " \
                  "SO=coordinate " \
                  "RGID=1 " \
                  "RGLB=Lib1 " \
                  "RGPL=illumina " \
                  "RGPU=hiseq2000 " \
                  "RGSM=%s" % (input_file,
                               output_file,
                               project_name)
    return cmd_replace


def remove_duplicates(inbamfile,
                      outbamfile,
                      metrics_file):
    """
    Remove duplicate reads.
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file
    :param metrics_file: file with summary statistics about
    :return: Command to be executed; type str
    """

    cmd_rmdup = "java -Xmx6g -jar $NGS_PICARD/MarkDuplicates.jar " \
                "INPUT=%s " \
                "OUTPUT=%s " \
                "METRICS_FILE=%s " \
                "REMOVE_DUPLICATES=true" % (inbamfile,
                                            outbamfile,
                                            metrics_file)
    return cmd_rmdup


def index_bam(inbamfile):
    """
    Index BAM file.
    :param inbamfile: name of BAM formatted file
    :return: Command to be executed; type str
    """

    cmd_index_bam = "samtools index %s" % (inbamfile)
    return cmd_index_bam


def splitntrim(input_file,
               output_file,
               ref_genome):
    """
    source: http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq
    - SplitNCigarReads developed specially for RNAseq,
    which splits reads into exon segments (getting
    rid of Ns but maintaining grouping information)
    and hard-clip any sequences overhanging into the
    intronic regions. -
    :param project_name: name of project (given by user
    :param output_dir: where the output files should be written
    :param ref_genome: reference genome  (.fa)
    :return: message to be logged & command to be executed; type str
    """

    cmd_splitntrim = "java -Xmx6g -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                     "-T SplitNCigarReads " \
                     "-R %s.fa " \
                     "-I %s " \
                     "-o %s " \
                     "-rf ReassignOneMappingQuality " \
                     "-RMQF 255 " \
                     "-RMQT 60 " \
                     "-U ALLOW_N_CIGAR_READS" % (ref_genome,
                                                 input_file,
                                                 output_file)
    return cmd_splitntrim


def bqsr(input_file,
         output_file,
         ref_genome,
         recal_report):
    """
    source: http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq
    - We do recommend running base recalibration (BQSR). Even though the
    effect is also marginal when applied to good quality data, it can 
    absolutely save your butt in cases where the qualities have 
    systematic error modes.
    Both steps 4 and 5 are run as described for DNAseq (with the same 
    known sites resource files), without any special arguments. Finally, 
    please note that you should NOT run ReduceReads on your RNAseq data. 
    The ReduceReads tool will no longer be available in GATK 3.0. -
    """

    cmd_splitntrim = "java -Xmx6g -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                     "-T PrintReads " \
                     "-R %s.fa " \
                     "-I %s " \
                     "-BQSR %s.recal_data.grp " \
                     "-o %s " % (ref_genome,
                                 input_file,
                                 recal_report,
                                 output_file)
    return cmd_splitntrim


def varcall_bamfo(input_file,
                  output_file_gatk,
                  ref_genome):
    """
    Variant calling with bamfo (author: Tomas Konopka).
    """
    cmd_bamfo = "bamfo callvariants " + " \\\n" \
                + "--bam " + input_file + " \\\n" \
                + "--output " + output_file_gatk + " \\\n" \
                + "--genome " + ref_genome

    return cmd_bamfo


def varcall_samtools(input_file,
                     output_file,
                     ref_genome):
    """
    Variant calling with samtools.
    """

    cmd_samtools = "samtools mpileup " \
                   "-C 50 " \
                   "-uf %s " \
                   "%s " \
                   "| bcftools view -vcg - > %s" % (ref_genome,
                                                    input_file,
                                                    output_file)
    return cmd_samtools


def varcall_gatk(input_file,
                 output_file_gatk,
                 ref_genome):
    """
    http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq
    'We have added some functionality to the variant calling code which
    will intelligently take into account the information about intron-exon
    split regions that is embedded in the BAM file by SplitNCigarReads.
    In brief, the new code will perform dangling head merging operations
    and avoid using soft-clipped bases (this is a temporary solution) as
    necessary to minimize false positive and false negative calls. To invoke
    this new functionality, just add -dontUseSoftClippedBases to your regular
    HC command line. Note that the -recoverDanglingHeads argument which was
    previously required is no longer necessary as that behavior is now enabled
    by default in HaplotypeCaller. Also, we found that we get better results
    if we lower the minimum phred-scaled confidence threshold for calling
    variants on RNAseq data, so we use a default of 20 (instead of 30 in
    DNA-seq data).'
    """

    cmd_varcall_gatk = "java -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                       "-T HaplotypeCaller " \
                       "-R %s.fa " \
                       "-I %s " \
                       "-dontUseSoftClippedBases " \
                       "-stand_call_conf 20.0 " \
                       "-stand_emit_conf 20.0 " \
                       "-o %s" % (ref_genome,
                                  input_file,
                                  output_file_gatk)
    return cmd_varcall_gatk


def variant_filtering(input_file,
                      output_file,
                      ref_genome):
    """
    To filter the resulting callset, you will need to apply hard filters,
    as we do not yet have the RNAseq training/truth resources that would
    be needed to run variant recalibration (VQSR).
    We recommend that you filter clusters of at least 3 SNPs that are within
    a window of 35 bases between them by adding -window 35 -cluster 3 to your
    command. This filter recommendation is specific for RNA-seq data.
    As in DNA-seq, we recommend filtering based on Fisher Strand values
    (FS > 30.0) and Qual By Depth values (QD < 2.0).
    """

    cmd_filter = "java -jar $NGS_GATK/GenomeAnalysisTK.jar " \
                 "-T VariantFiltration " \
                 "-R %s.fa " \
                 "-V %s " \
                 "-window 35 " \
                 "-cluster 3 " \
                 "-filterName FS " \
                 "-filter \"FS > 30.0\" " \
                 "-filterName QD " \
                 "-filter \"QD < 2.0\" " \
                 "-o %s" % (ref_genome,
                            input_file,
                            output_file)
    return cmd_filter


def convert2annovar(input_file,
                    output_file):
    """
    Annotate vcf file with annovar.
    """

    cmd_annovar = "convert2annovar.pl " \
                  "--format vcf4 " \
                  "--includeinfo %s >" \
                  "%s.annovar" % (input_file,
                                 output_file)
    return cmd_annovar


def dbsnp_filter(dbtype,
                 buildversion,
                 input_file,
                 annovar_dir):

    cmd_dbsnp_filter = "annotate_variation.pl " \
                        "--filter "\
                        "--dbtype %s " \
                        "--buildver %s " \
                        "%s %s" % (dbtype,
                                   buildversion,
                                   input_file,
                                   annovar_dir)
    return cmd_dbsnp_filter


def gene_annotation(buildversion,
                    input_file,
                    annovar_dir):

    cmd_dbsnp_filter = "annotate_variation.pl " \
                        "--buildver %s " \
                        "%s %s" % (buildversion,
                                   input_file,
                                   annovar_dir)
    return cmd_dbsnp_filter
