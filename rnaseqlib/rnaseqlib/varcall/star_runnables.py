"""
Collection of functions that execute external software or UNIX commands.
"""


def rnaseq_align(star_genome,
                 gtf,
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
                "--sjdbGTFfile %s " \
                "--readFilesIn %s %s " \
                "--runThreadN %s " \
                "--genomeLoad NoSharedMemory " \
                "--outFilterIntronMotifs RemoveNoncanonical " \
                "--outSAMtype BAM SortedByCoordinate " \
                "--outFileNamePrefix %s" % (star_genome, gtf, read1,
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
    --sjdbOverhang
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
