from __future__ import division
from Bio import pairwise2


def _jfill(j_ref, seq):
    """ Replace any nucleotides missing from the start of the J gene
    """

    j_ref = j_ref[-25:]
    j_seq = seq[-25:]

    j_filled = _fill_missing_nt(j_ref, j_seq)

    seq_filled = seq[:-25] + j_filled

    return seq_filled


def _vfill(v_ref, seq):
    """ Replace any nucleotides missing from the start of the V gene
    """

    v_ref = v_ref[:75]
    v_seq = seq[:75]

    v_filled = _fill_missing_nt(v_ref, v_seq)
    seq_filled = v_filled + seq[75:]

    return seq_filled


def global_pw_align(s0, s1):
    """ Wrapper for biopython global pairwise alignment.
        Has a strong gap penalty except end gaps are not penalized.
    """

    aln = pairwise2.align.globalms(s0, s1, 2, -1, -6, -1,
                                   penalize_end_gaps=False)[0]
    s0_aln = aln[0]
    s1_aln = aln[1]

    assert len(s0_aln) == len(s1_aln)

    return s0_aln, s1_aln


def _fill_missing_nt(ref_seq, seq):
    """ Perform pairwise alignment between two sequences and fill in
        missing nucleotides
    """

    ref_aligned, seq_aligned = global_pw_align(ref_seq, seq)

    if '-' in seq_aligned:
        filled = ''
        for char_ref, char_seq in zip(ref_aligned, seq_aligned):
            if char_seq == '-':
                filled += char_ref
            else:
                filled += char_seq
    else:
        filled = seq_aligned

    return filled


def fill_clean_sequence(seq_id, v_ref_seq, j_ref_seq, seq):
    """ Clean gaps arising from deletions and fill in missing nucleotides at
        the beginning of the V gene or end of the J gene.
    """

    if '-' in seq:
        degap = seq.replace('-','')
    else:
        degap = seq

    vfilled = _vfill(v_ref_seq, degap)

    vjfilled = _jfill(j_ref_seq, vfilled)

    if seq != vjfilled:
        s0, s1 = global_pw_align(seq, vjfilled)
        s0 = '{}...{}'.format(s0[:25], s0[-25:])
        s1 = '{}...{}'.format(s1[:25], s1[-25:])
        print('{}:\n\tBefore:\t{}\n\tFilled:\t{}'.format(seq_id, s0, s1))

    return vjfilled
