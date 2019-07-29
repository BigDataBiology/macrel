#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastformatterCommandline

file1 = sys.argv[1]
prepbed = sys.argv[2]
leng = sys.argv[3]

def create_and_run_bedrequest(blast_xml_filename, bed_out_filename, qlen):
    '''
    Take a tblastn output (outputfmt 7), identify ORFS of interest and
    create a file 'bedrequests.bed', which can be used by bedtools to
    return the actual nucleotide sequences of our matches.

    Parameters
    ----------
    blast_xml_filename : str
    bed_out_filename : str
    qlen : int
    '''

    # tblastn returns an alignment of the translated dna sequences it found. We, however,
    # want that alignment PLUS all the dna that is around it. (that is, everything that is
    # black in the pictures.) So, we take every result, check where it came from (e.g. position
    # 123456 - 123600 in genome XYZ.). Then, we extract a dna that comes right before and after
    # that, e.g. (123400 - 123465) and (123591 - 123640). Depending on how long the homologue
    # protein is compared to the query protein, these "appendices" are shorter or longer.

    # In mayor step 3, the appendices and the homologue proteins are then glued together to
    # receive our pictures.

    # We want appendices to overlap with our homologue protein, so we can check
    # everything is right.
    secoverlap = 3

    healthy_seqs = {}

    # Load xml files created by tblastn. The xml files contain the information
    # about where, who and what the homologue proteins are.
    NcbiblastformatterCommandline(archive=blast_xml_filename,
                                  outfmt='7 qseqid sseqid evalue qseq sseq',
                                  out="{}.txt".format(blast_xml_filename.split('.')[0]))
    # For every matched sequence...
    with open(blast_xml_filename) as blast_file:
        blast_records = NCBIXML.parse(blast_file)
        with open(bed_out_filename, "w") as bed_file:
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        # Update list of known sequences
                        healthy_seqs[alignment.title] = hsp.sbjct

                        # Extract the critical information needed for bedrequest

                        # a) where in the genome is the homologue?
                        nuc_start = hsp.sbjct_start
                        nuc_end = hsp.sbjct_end
                        
                        # Plus or minus strand?
                        strand = '+' if (hsp.frame[1] > 0) else '-'

                        # b) which parts of the query does this homologue cover?
                        # (e.g. amino acids 9 to 104)
                        q_start = hsp.query_start
                        q_end = hsp.query_end

                        startorig = min(nuc_start, nuc_end)
                        stoporig = max(nuc_start, nuc_end)

                        # How far do we have to reach back and forth, to cover the
                        # whole query protein (plus 20aa at the end?)
                        walkback = int(q_start) - 1
                        walkfront = int(qlen) - (q_end)
			             
                        # Also, we want our 3aa overlap.
                        extend = 1
                        if extend == 1:
                            if strand == '+':
                                start = startorig - (3 * walkback)
                                stop = stoporig + (3 * walkfront)
                            elif strand == '-':
                                start = startorig - (3 * walkfront)
                                stop = stoporig + (3 * walkback)

                        # and write ID, Start, Stop and Strand (+/-) to
                        # bedrequests.bed, which will later be used by
                        # 'bedtools getfasta ...'

                        bed_out_template = '{}\t{{}}\t{{}}\t{{}}\t{{}}\t{{}}'.format(alignment.hit_id)
                        bed_out_lines = []
                        vals = [(start - 1, startorig - 1 + (3 * secoverlap)),
                                (stoporig - (3 * secoverlap), stop)]
                        tail = ["{}|{}:{}|{}|{}".format(alignment.hit_id,
                                                        q_start,
                                                        q_end,
                                                        hsp.sbjct,
                                                        hsp.expect),
                                'dummy',
                                strand]
                        for val1, val2 in vals:
                            line = bed_out_template.format(val1,
                                                           val2, *tail)
                            bed_out_lines.append(line)
                        if strand == '-':
                            bed_out_lines.reverse()
                        
                        if (vals[0][0] > 0) and (vals[0][1] > 0):
                            for line in bed_out_lines:
                                print(line, file=bed_file)

create_and_run_bedrequest(file1, prepbed, leng)
