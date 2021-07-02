
"""
a python version of sam2tsv from jarkit (http://lindenb.github.io/jvarkit/Sam2Tsv.html)

for testing
sam_file = pysam.AlignmentFile("/home/zh/local/src/samtools-1.3/examples/toy.sam", "r", check_sq=False)
ref_fa = pysam.FastaFile("/home/zh/local/src/samtools-1.3/examples/toy.fa")

=================================
source of code: https://gist.github.com/mt1022/737ef20f43d5acd4bc75dba0be8f334b
author: mt1022 (zh) altered by Leanne S Whitmore
date: 2017-11-28 22:39
"""
import sys
import itertools
import pysam

def sam2tsv_function(alignfile, ref_file, variant_pos, outputfile="tmp.txt"):
    '''modified code from mt10022 to generate this function'''

    cigar_dict = {i:j for i, j in zip(range(10), 'MIDNSHP=XB')}
    op_consume_query = [0, 1, 4, 7, 8]
    variant_pos = variant_pos -1
    sam_file = pysam.AlignmentFile(alignfile, 'rb')
    if alignfile == 'sam':
        sam_file = pysam.AlignmentFile(alignfile, 'r')

    ref_fa = pysam.FastaFile(ref_file)
    outputf = open(outputfile, "w")
    for align in sam_file:
        #  if align.query_name in reads:
        ref_name = sam_file.getrname(align.reference_id)
        query_seq = align.query_sequence
        query_qualities = align.query_qualities
        ref_pos = align.get_reference_positions(full_length=True)

        # only include cigars that consume query
        cigar_lolist = [[cigar_dict[i]]*j for i, j in align.cigartuples if i in op_consume_query]
        cigars = list(itertools.chain.from_iterable(cigar_lolist))

        read_info = [align.query_name, align.flag, ref_name]
        for i in range(align.query_length):
            query_qual = '.' if query_qualities is None else query_qualities[i]
            # None values will be included for any soft-clipped or unaligned positions within the read.
            if ref_pos[i] is None:
                out = read_info + [i, query_seq[i], query_qual, '.', '.', cigars[i]]
            else:
                try:
                    ref_base = ref_fa.fetch(reference=ref_name, start=ref_pos[i], end=ref_pos[i]+1)
                except KeyError:
                    ref_base = '-'
                out = read_info + [i, query_seq[i], query_qual, ref_pos[i] + 1, ref_base, cigars[i]]
            if ref_pos[i]==variant_pos:
                outputf.write('\t'.join([str(j) for j in out])+"\n")
    outputf.close()