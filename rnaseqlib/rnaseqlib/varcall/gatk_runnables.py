"""
Collection of functions that execute external software or UNIX commands.
"""


def splitntrim(input_file,
               output_file,
               ref_genome,
               heap_mem):
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

    cmd_splitntrim = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                     "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
                     "-T SplitNCigarReads " \
                     "-R %s.fa " \
                     "-I %s " \
                     "-o %s " \
                     "-rf ReassignOneMappingQuality " \
                     "-RMQF 255 " \
                     "-RMQT 60 " \
                     "-U ALLOW_N_CIGAR_READS" % (heap_mem,
                                                 ref_genome,
                                                 input_file,
                                                 output_file)
    return cmd_splitntrim


def target_creator(input_file,
                   output_file,
                   ref_genome,
                   known_1000g,
                   known_mills,
                   heap_mem):
    """
    Perform local realignment around indels to correct mapping-related artifacts.
    """

    cmd_target = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                 "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
                 "-T RealignerTargetCreator " \
                 "-R %s.fa " \
                 "-I %s " \
                 "-known %s " \
                 "-known %s " \
                 "-o %s" % (heap_mem,
                            ref_genome,
                            input_file,
                            known_1000g,
                            known_mills,
                            output_file)
    return cmd_target


def indel_realignment(input_file,
                      output_file,
                      ref_genome,
                      known_1000g,
                      known_mills,
                      target,
                      heap_mem):


    cmd_indel = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
                "-T IndelRealigner " \
                "-R %s.fa " \
                "-I %s " \
                "-targetIntervals %s " \
                "-known %s " \
                "-known %s " \
                "-o %s" % (heap_mem,
                           ref_genome,
                           input_file,
                           target,
                           known_1000g,
                           known_mills,
                           output_file)

    return cmd_indel


def analyse_covariation_patterns(input_file,
                                 output_file,
                                 ref_genome,
                                 known_1000g,
                                 known_mills,
                                 dbsnp,
                                 heap_mem):
    """
    Analyze patterns of covariation in the sequence dataset.
    """

    cmd_acp = "java -%s -Djava.io.tmpdir=$TMPDIR " \
              "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
              "-T BaseRecalibrator " \
              "-R %s.fa " \
              "-I %s " \
              "-knownSites %s " \
              "-knownSites %s " \
              "-knownSites %s " \
              "-o %s" % (heap_mem,
                         ref_genome,
                         input_file,
                         known_1000g,
                         known_mills,
                         dbsnp,
                         output_file)
    return cmd_acp


def analyse_covariation_patterns_2ndpass(input_file,
                                         output_file,
                                         ref_genome,
                                         known_1000g,
                                         known_mills,
                                         dbsnp,
                                         recal,
                                         heap_mem):
    """
    Do a second pass to analyze covariation remaining after recalibration.
    """

    cmd_acp2ndpass = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                     "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
                     "-T BaseRecalibrator " \
                     "-R %s.fa " \
                     "-I %s " \
                     "-knownSites %s " \
                     "-knownSites %s " \
                     "-knownSites %s " \
                     "-BQSR %s " \
                     "-o %s" % (heap_mem,
                                ref_genome,
                                input_file,
                                known_1000g,
                                known_mills,
                                dbsnp,
                                recal,
                                output_file)
    return cmd_acp2ndpass


def plot_recalibration(ref_genome,
                       before,
                       after,
                       plot_name,
                       heap_mem):
    """
    Generate before/after plots
    """

    cmd_plot = "java -%s -Djava.io.tmpdir=$TMPDIR " \
               "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
               "-T AnalyzeCovariates " \
               "-R %s.fa " \
               "-before %s " \
               "-after %s " \
               "-plots %s " % (heap_mem,
                               ref_genome,
                               before,
                               after,
                               plot_name)
    return cmd_plot


def bqsr(ref_genome,
         input_file,
         bqsr,
         output_file,
         heap_mem):

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

    cmd_bqsr = "java -%s -Djava.io.tmpdir=$TMPDIR " \
               "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
               "-T PrintReads " \
               "-R %s.fa " \
               "-I %s " \
               "-BQSR %s " \
               "-o %s " % (heap_mem,
                           ref_genome,
                           input_file,
                           bqsr,
                           output_file)
    return cmd_bqsr


def varcall_gatk(input_file,
                 output_file_gatk,
                 ref_genome,
                 heap_mem):
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

    cmd_varcall_gatk = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                       "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
                       "-T HaplotypeCaller " \
                       "-R %s.fa " \
                       "-I %s " \
                       "-dontUseSoftClippedBases " \
                       "-stand_call_conf 20.0 " \
                       "-stand_emit_conf 20.0 " \
                       "-o %s" % (heap_mem,
                                  ref_genome,
                                  input_file,
                                  output_file_gatk)
    return cmd_varcall_gatk


def variant_filtering(input_file,
                      output_file,
                      ref_genome,
                      heap_mem):
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

    cmd_filter = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                 "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
                 "-T VariantFiltration " \
                 "-R %s.fa " \
                 "-V %s " \
                 "-window 35 " \
                 "-cluster 3 " \
                 "-filterName FS " \
                 "-filter \"FS > 30.0\" " \
                 "-filterName QD " \
                 "-filter \"QD < 2.0\" " \
                 "-o %s" % (heap_mem,
                            ref_genome,
                            input_file,
                            output_file)
    return cmd_filter


def var2table_rna(input_file,
                   output_file,
                   ref_genome,
                   heap_mem):
    """
    Extract fields from VCF file.
    """

    cmd_2table = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                 "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
                 "-R %s " \
                 "-T VariantsToTable " \
                 "-V %s " \
                 "-F CHROM " \
                 "-F POS " \
                 "-F REF " \
                 "-F ALT " \
                 "-F QUAL " \
                 "-F DP " \
                 "-F QD " \
                 "-F Func.ensGene " \
                 "-F Gene.ensGene " \
                 "-F ExonicFunc.ensGene " \
                 "-F AAChange.ensGene " \
                 "-GF GT " \
                 "-GF AD " \
                 "-F SIFT_score " \
                 "-F Polyphen2_HDIV_score " \
                 "-F CADD_phred " \
                 "-o %s " \
                 "--allowMissingData" % (heap_mem,
                                         ref_genome,
                                         input_file,
                                         output_file)
    return cmd_2table


def var2table_dna(input_file,
                  output_file,
                  ref_genome,
                  heap_mem):
    """
    Extract fields from VCF file for dna (exomes in particular; in addition
    these also include vqslod score and SNPeff annotation
    """

    cmd_2table = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                 "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
                 "-R %s " \
                 "-T VariantsToTable " \
                 "-V %s " \
                 "-F CHROM " \
                 "-F POS " \
                 "-F ID " \
                 "-F REF " \
                 "-F ALT " \
                 "-F SNPEFF_GENE_NAME " \
                 "-F Gene.ensGene " \
                 "-F SNPEFF_GENE_BIOTYPE " \
                 "-F Func.ensGene " \
                 "-F ExonicFunc.ensGene " \
                 "-F AAChange.ensGene " \
                 "-F QUAL " \
                 "-F QD " \
                 "-GF GT " \
                 "-GF AD " \
                 "-GF DP " \
                 "-GF GQ " \
                 "-GF PL " \
                 "-F VQSLOD " \
                 "-F cytoBand " \
                 "-F genomicSuperDups " \
                 "-F snp129 " \
                 "-F snp142Common " \
                 "-F 1000g2015feb_all " \
                 "-F 1000g2015feb_eur " \
                 "-F 1000g2015feb_afr " \
                 "-F 1000g2015feb_amr " \
                 "-F 1000g2015feb_eas " \
                 "-F 1000g2015feb_sas " \
                 "-F esp5400_all " \
                 "-F esp6500siv2_all " \
                 "-F cosmic70 " \
                 "-F clinvar_20150330 " \
                 "-F SIFT_score " \
                 "-F SIFT_pred " \
                 "-F Polyphen2_HDIV_score " \
                 "-F Polyphen2_HDIV_pred " \
                 "-F Polyphen2_HVAR_score " \
                 "-F Polyphen2_HVAR_pred " \
                 "-F LRT_score " \
                 "-F LRT_pred " \
                 "-F MutationTaster_score " \
                 "-F MutationTaster_pred " \
                 "-F FATHMM_score " \
                 "-F FATHMM_pred " \
                 "-F RadialSVM_score " \
                 "-F RadialSVM_pred " \
                 "-F LR_score " \
                 "-F LR_pred " \
                 "-F VEST3_score " \
                 "-F CADD_raw " \
                 "-F CADD_phred " \
                 "-F GERP++_RS " \
                 "-F phyloP46way_placental " \
                 "-F phyloP100way_vertebrate " \
                 "-F SiPhy_29way_logOdds " \
                 "-F caddgt10 " \
                 "-F CADDindel " \
                 "-F CADDindel_Phred " \
                 "-o %s " \
                 "--allowMissingData" % (heap_mem,
                                         ref_genome,
                                         input_file,
                                         output_file)
    return cmd_2table


def var2table_proudpv(input_file,
                      output_file,
                      ref_genome,
                      heap_mem):
    """
    Extract fields from VCF file.
    """

    cmd_2table = "java -%s -Djava.io.tmpdir=$TMPDIR " \
                 "-jar $NGS_GATK/GenomeAnalysisTK.jar " \
                 "-R %s " \
                 "-T VariantsToTable " \
                 "-V %s " \
                 "-F CSQT " \
                 "-F Func.ensGene " \
                 "-F ExonicFunc.ensGene " \
                 "-F CHROM " \
                 "-F POS " \
                 "-F REF " \
                 "-F ALT " \
                 "-F AAChange.ensGene " \
                 "-F DP " \
                 "-GF GQ " \
                 "-GF AD " \
                 "-GF VF " \
                 "-F Gene.ensGene " \
                 "-F cosmic70 " \
                 "-F SIFT_score " \
                 "-F CADD_phred " \
                 "-F Polyphen2_HDIV_score " \
                 "-o %s " \
                 "--allowMissingData" % (heap_mem,
                                         ref_genome,
                                         input_file,
                                         output_file)
    return cmd_2table