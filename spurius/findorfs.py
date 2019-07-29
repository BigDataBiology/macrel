'''This file contains most of the important functions and things.
It is imported by executer.py, and the functions are executed
throughout the run of executer.py.'''

import argparse
import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastformatterCommandline

access = sys.argv[1]
input = sys.argv[2]
leng = sys.argv[3]
nmb = sys.argv[4]
qt = sys.argv[5]
ev = sys.argv[6]
DICIONARIO = sys.argv[7]
outdir = 'npdir'


TBLASTN_COVERED = 4
TBLASTN_GAP = 0
TBLASTN_START_CODON = 99
ANYWHERE_STOP_CODON = 300
FRONT_STOP = 1
BODY_STOP = 2
BACK_STOP = 3    
ALIGNMENT_BODY = 1
PSEUDOCOUNT = 1

def load_sequence_files(fasta_filename):
    """Load FASTA file.

    If the FASTA file contains an odd number of sequences for a particular
    identifier, they are all discarded.

    Parameters
    ----------
    str
        The filename of a FASTA file.

    Returns
    -------
    [Bio.SeqRecord.SeqRecord]
    """
    files = list(SeqIO.parse(fasta_filename, 'fasta'))
    if len(files) < 2:
        raise RuntimeError('Only 0 or 1 tblastn hits found.')

    # Changing name. This became necessary after switching from bedtools version 2.20 to 2.26
    for file_ in files:
        file_.id = file_.id.split('::')[0]

    # Remove single entries. They stem from tblastn results too close to genome borders.

    last_filename = None
    first_n_of_filename = 0
    n_file = 0
    bad_ranges = []
    for i, file_ in enumerate(files):
        filename = file_.id
        if filename == last_filename:
            n_file += 1
        else:
            if n_file % 2 != 0:
                print('Warning 0: Found an issue in one of the two '
                      ' appendices to sequence {}. Removing {} '
                      'sequences.'.format(i, first_n_of_filename - i))
                bad_ranges.append((first_n_of_filename, i))
            first_n_of_filename = i
            n_file = 1
        last_filename = filename
    bad_ranges.reverse()
    for start, end in bad_ranges:
        del files[start:end]

    # A: append our bedtool results to tblast output. This is a sensible step.
    # Before merging, we make several quality checks to assure that the
    # sequences look exactly as they shoud.

    assert len(files) % 2 == 0, ('Odd number of sequences. At least one lonely '
                                 'sequence was not captured.')
    return files

def combine_sequences(sequences):
    """Combine pairs of sequences found by TBLASTN.
    Parameters
    ----------
    [Bio.SeqRecord.SeqRecord]

    Returns
    -------
    [str]
    Raises
    ------
    RuntimeError
        If the sequences in the file can't be combined.
    """

    newseq = []
    for i in range(0, len(sequences), 2):
        file_ = sequences[i]
        next_file = sequences[i + 1]
        if len(file_) % 3:
            raise RuntimeError('Length of front appendix '
                               'is not a multiple of 3')

        if len(next_file) % 3:
            raise RuntimeError('Length of end appendix is not a multiple of 3')

        if file_.id != next_file.id:
            raise RuntimeError('Front and End appendix do '
                               'not have the same ID.')

        front = str(file_.seq.translate())
        body = file_.id.split('|')[-2]
        back = str(next_file.seq.translate())
        front_tail = front[-3:]
        body_head = body[:3]
        body_tail = body[-3:]
        back_head = back[:3]

        if not (_seq_contains_gaps_or_pseudocode(front_tail)
                or _seq_contains_gaps_or_pseudocode(body_head)
                or front_tail == body_head):
            raise RuntimeError('Seq {}: Front appendix and MSA body do not'
                               ' overlap correctly. (what should be there {} '
                               'vs what tblastn claims {}), aka {} <- what is'
                               ' really in the genome(bedtools)'.format(
                                   i, front_tail, body_head, file_.seq[-9:]))
        if not (_seq_contains_gaps_or_pseudocode(body_tail)
                or _seq_contains_gaps_or_pseudocode(back_head)
                or body_tail == back_head):
            raise RuntimeError('Seq {}: End appendix and MSA body do not '
                               'overlap correctly. (what should be there {} '
                               'vs what tblastn claims {}), aka {}  <- what '
                               'is really in the genome(bedtools)'.format(
                                   i, back_head, body_tail, next_file.seq[0:9]))


        # If everything seems fine, merge! This is the important step #
        # Just to be clear, we are taking the homologue sequence found by
        # tblastn as-is, including gaps.
        # We then glue our appendices to it. The reason we do it so strangely
        # is that we need the gaps.
        newseq.append((file_.id, front[:-3], body, back[3:]))
    return newseq

def make_visualisation_matrix(combined_sequences, qlen):
    n_seq = len(combined_sequences)
    t = int(qlen) + 20
    pic = np.zeros((n_seq, t))
    pic_stops = np.zeros((n_seq, t))
    counter = 0
    all_stop_codons = []
    sseqlens = np.zeros(n_seq)
    starts = np.zeros(n_seq)
    stops = np.zeros(n_seq)
    evals = np.zeros(n_seq)
    for file_id, front_tail, body, back_head in combined_sequences:
        *_, protein_coords, seq, evalue_s = file_id.split("|")
        # Fix formatting of evalues from bedtools in some circumstances
        evalue_s = evalue_s.split('(')[0]
        start, stop = map(int, protein_coords.split(':'))
        evals[counter] = float(evalue_s)
        if evals[counter] == 0:
            evals[counter] = 1e-200
        starts[counter] = start
        stops[counter] = stop
        pic[counter, start - 1:stop] = -np.log10(evals[counter])
        pic_stops[counter, start - 1:stop] = TBLASTN_COVERED
        # Appendix before seq
        front_stop_codons = [i for i, letter in enumerate(front_tail)
                             if letter == '*']
        front_start_codons = [i for i, letter in enumerate(front_tail)
                              if letter == 'M']

        # Appendix after seq
        back_stop_codons = [i + stop - 1 for i, letter in enumerate(back_head)
                            if letter == '*']
        back_start_codons = [i + stop - 1 for i, letter in enumerate(back_head)
                             if letter == 'M']

        # Actual seq
        sseqlen = len(seq)
        sseqlen_no_gaps = len(seq.replace('-', ''))
        qseqlen = (stop - start) + 1

        body_stop_codons = [round((j / sseqlen) * qseqlen) + start
                            for j, letter in enumerate(body) if letter == '*']
        body_start_codons = [round((j / sseqlen) * qseqlen) + start
                             for j, letter in enumerate(body) if letter == 'M']
        body_gaps = [round((j / sseqlen) * qseqlen) + start
                     for j, letter in enumerate(body) if letter == '-']

        all_stop_codons = front_stop_codons + back_stop_codons + body_stop_codons

        # For our histograms, we only accept "post-alignment stops" if the
        # alignment before covered at least 75% of the query.
        #if (len(files[i+1].seq.translate()[3:])-20) < (0.25 * (len(file_.id.split('|',5)[-2])-20)):
        #    occs_late.append(back_stop_codons)

        pic[counter, body_gaps] = TBLASTN_GAP
        pic[counter, body_start_codons] = TBLASTN_START_CODON
        pic[counter, front_start_codons] = TBLASTN_START_CODON
        pic[counter, back_start_codons] = TBLASTN_START_CODON
        pic[counter, all_stop_codons] = ANYWHERE_STOP_CODON
        pic_stops[counter, body_stop_codons] = BODY_STOP
        pic_stops[counter, back_stop_codons] = BACK_STOP
        pic_stops[counter, front_stop_codons] = FRONT_STOP
        sseqlens[counter] = sseqlen_no_gaps
        counter += 1
        all_stop_codons = []

    ###################################################################################################################################
    ### Now that we have our merged sequences and identified stop codons and evalues, we do a bit of shuffling and annoying stuff.  ###
    ### It might even be that we dont need some of these lines anymore. I can look through it when I have time.                     ###
    ###################################################################################################################################

    # Reorder sequences by evalue
    evals = np.array(evals)
    evals_sort = np.argsort(-1 * evals)
    sorted_pic = pic[evals_sort]#[np.argsort(np.max(pic2, axis = i))]
    sorted_pic_stops = pic_stops[evals_sort]
    starts = starts[evals_sort]
    stops = stops[evals_sort]
    sseqlens = sseqlens[evals_sort]

    return evals, evals_sort, sorted_pic, sorted_pic_stops, starts, stops, sseqlens


def count_stops_per_aa(mat_stops, starts, stops):
    '''
    Input: The colorful MSA in matrix form
    Output: Number of Stop codons per amino acid.
    '''
    # Remove 20 aa appended to end
    mat_stops = mat_stops[:,:-20]
 
    # Number of amino acids to ignore at each start and end
    aa_tolerance = 10   
    
    # We don't even need to continue, if there are no homologous seqs or the 
    # sequence is smaller than 2 times the tolerance. 
    if ((np.shape(mat_stops)[0] == 0) or 
                      (np.shape(mat_stops)[1] <= 2*aa_tolerance)):
        return 0 + PSEUDOCOUNT, 0


    
    # Body_mask is an array of 0 and 1, where 1 marks amino acids in the body.
    body_mask = np.zeros(np.shape(mat_stops))
    for hom_seq in range(len(mat_stops)):
        body_mask[hom_seq, int(starts[hom_seq] + aa_tolerance - 1)
                         : int(stops[hom_seq] - aa_tolerance)] = ALIGNMENT_BODY
    
    mat_masked = mat_stops * body_mask
    
    n_stop = np.sum(mat_masked == BODY_STOP)
    n_bodyaas = np.sum(mat_masked > 0) + PSEUDOCOUNT
    stops_per_aa = n_stop / n_bodyaas
    
    return n_bodyaas, stops_per_aa


def find_stops(inputfile, qname, n_fasta, qlen, querytype, eval_cutoff, np_dir):
        # Load nucleotide appendices. The MSA body is in the ID.
    
    files = load_sequence_files(inputfile)
    newseq = combine_sequences(files)
    evals, evals_sort, pic2, pic_stops, starts, stops, sseqlens = make_visualisation_matrix(newseq, qlen)
    stop_masked_pic = np.ma.masked_where(pic2 != ANYWHERE_STOP_CODON, pic2)

    # Mask out elements which don't contain stop codons and then set them to be grey
    start_masked_pic = np.ma.masked_where(pic2 != TBLASTN_START_CODON, pic2)
    start_masked_pic[np.ma.nonzero(start_masked_pic)] = 40

    pic2_onlycol = np.ma.masked_where((pic2 == TBLASTN_GAP) |
                                      (pic2 == ANYWHERE_STOP_CODON) |
                                      (pic2 == TBLASTN_START_CODON), pic2)


    # Get the "Spurious_ORF_XY" name right...
    if querytype == "af":
        try:
            a = np.load(DICIONARIO)
            truename = qname.split(' ')[0].split('|')[-1]
            spuriosity = a[np.where(a[:, 0] == truename)[0][0],1]
        except:
            spuriosity = 'unknown_af'
    else:
        spuriosity = 'noafhit'

    calc = -np.log10(evals[evals_sort])
    idxs_above_th = [ calc > int(eval_cutoff) ]
    pic_stops_cutoff = pic_stops[idxs_above_th]

    sseqlens = sseqlens[idxs_above_th]
    starts = starts[idxs_above_th]
    stops = stops[idxs_above_th]
    sseqlens = np.reshape(sseqlens, (len(sseqlens), 1))
    starts = np.reshape(starts, (len(starts), 1))
    stops = np.reshape(stops, (len(stops), 1))
    n_homologous_seqs = np.shape(pic_stops_cutoff)[0]
    if np.shape(pic_stops_cutoff)[0] == 0:
        raise RuntimeError('No sequences left after '
                           'evalue cutoff. Earlier: {}'.format(
                               np.shape(pic_stops)[0]))

    # Extract the number of stop codons in the region of interest
    body_aa, stops_per_aa = count_stops_per_aa(pic_stops_cutoff, starts, stops)

    #Save the colorful MSA
    file_name_template = "{}_{}_{}_{}.{{}}".format(querytype,
                                                   n_fasta,
                                                   spuriosity,
                                                   qname)

    # Save a much reduced picture, containing only 0s and 1,2,3 at STOP
    # positions before (1), within (2) or after (3) the 'body' alignment.
    # This is used as input for simple_classifier.
    np.save(os.path.join(np_dir, file_name_template.format("npy")),
            np.hstack((stops_per_aa, n_homologous_seqs, qlen, body_aa)))

    return np.array([qname.split('_')[-1]]), np.array([stops_per_aa, n_homologous_seqs, qlen, body_aa])

def _seq_contains_gaps_or_pseudocode(seq):
    return len(set(seq) & {'X', '-', '*'}) > 0

find_stops(input, access, nmb, leng, qt, ev, outdir)
