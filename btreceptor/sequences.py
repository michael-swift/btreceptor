from __future__ import division
from Bio import pairwise2
from Bio.Seq import Seq
from . import parsing
import os
import pandas as pd


# IMGT defined CDR lengths
# http://www.imgt.org/IMGTrepertoire/2D-3Dstruct/FClength/human/Hu_CDR_IGHKLV.html
# http://www.imgt.org/IMGTrepertoire/2D-3Dstruct/FClength/mouse/Mu_CDR_IGHKLV.html
cdr_12_aa_lengths = {
    'human': {
        'IGHV1': [(8, 8)],
        'IGHV2': [(10, 7)],
        'IGHV3': [(8, 7), (8, 8), (8, 10)],
        'IGHV4': [(8, 7), (9, 7), (10, 7)],
        'IGHV5': [(8, 8)],
        'IGHV6': [(10, 9)],
        'IGHV7': [(8, 8)],
        'IGKV1': [(6, 3)],
        'IGKV2': [(11, 3), (12, 3)],
        'IGKV3': [(6, 3), (7, 3)],
        'IGKV3D': [(6, 3), (7, 3)],
        'IGKV4': [(12, 3)],
        'IGKV5': [(6, 3)],
        'IGKV6': [(6, 3)],
        'IGLV1': [(8, 3), (9, 3)],
        'IGLV2': [(9, 3)],
        'IGLV3': [(6, 3)],
        'IGLV4': [(7, 7)],
        'IGLV5': [(9, 7)],
        'IGLV6': [(8, 3)],
        'IGLV7': [(9, 3)],
        'IGLV8': [(9, 3)],
        'IGLV9': [(7, 8)],
        'IGLV10': [(8, 3)]
    },
    'mouse': {
        'IGHV1': [(8, 8)],
        'IGHV2': [(8, 7)],
        'IGHV3': [(10, 7), (9, 7), (8, 7)],
        'IGHV4': [(8, 8)],
        'IGHV5': [(8, 8), (8, 8)],
        'IGHV6': [(8, 10)],
        'IGHV7': [(8, 10)],
        'IGHV8': [(10, 7)],
        'IGHV9': [(8, 8)],
        'IGHV10': [(8, 10)],
        'IGHV11': [(8, 8)],
        'IGHV12': [(10, 7), (9, 7)],
        'IGHV13': [(8, 10)],
        'IGHV14': [(8, 8)],
        'IGHV15': [(9, 8)],
        'IGHV16': [(9, 8)],
        'IGKV1': [(11, 3)],
        'IGKV2': [(11, 3)],
        'IGKV3': [(10, 3)],
        'IGKV4': [(7, 3), (5, 3)],
        'IGKV5': [(6, 3)],
        'IGKV6': [(6, 3)],
        'IGKV7': [(12, 3)],
        'IGKV8': [(12, 3)],
        'IGKV9': [(6, 3)],
        'IGKV10': [(6, 3)],
        'IGKV11': [(6, 3)],
        'IGKV12': [(6, 3)],
        'IGKV13': [(6, 3)],
        'IGKV14': [(6, 3)],
        'IGKV16': [(6, 3)],
        'IGKV17': [(6, 3)],
        'IGKV18': [(6, 3)],
        'IGKV19': [(6, 3)],
        'IGKV20': [(6, 3)],
        'IGLV1': [(9, 3)],
        'IGLV2': [(7, 7)]
    }
}


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


def fill_clean_sequence(row, vj_ref_dict, verbose):
    """ Clean gaps arising from deletions and fill in missing nucleotides at
        the beginning of the V gene or end of the J gene.

        Args:
            row (pd.Series): Series containg the following columns:
                - sequence_vdj
                - v_call
                - j_call
                - seqid
            vj_ref_dict (dict): key: v_call or j_call, values: sequence
                - example: {'IGHV1-18*01': 'CAGGTTC...'
    """

    if '-' in row.sequence_vdj:
        degap = row.sequence_vdj.replace('-', '')
    else:
        degap = row.sequence_vdj

    try:
        vfilled = _vfill(vj_ref_dict[row.v_call], degap)
    except KeyError:
        print('{} V gene was not found in the reference'.format(row.v_call))
        return None

    try:
        vjfilled = _jfill(vj_ref_dict[row.j_call], vfilled)
    except KeyError:
        print('{} J gene was not found in the reference'.format(row.j_call))
        return None

    if row.sequence_vdj != vjfilled and verbose:
        s0, s1 = global_pw_align(row.sequence_vdj, vjfilled)
        s0 = '{}...{}'.format(s0[:25], s0[-25:])
        s1 = '{}...{}'.format(s1[:25], s1[-25:])
        print('{}:\n\tBefore:\t{}\n\tFilled:\t{}'.format(row.seqid, s0, s1))

    return vjfilled


def translate(nt_seq):
    """ Returns translated nucleotide sequence """

    return str(Seq(nt_seq).translate())


def _no_stop_codon(aa_seq):
    """ Returns True if a sequence does not contain a stop codon,
        otherwise returns False
    """

    if '*' not in aa_seq:
        return True
    else:
        return False


def _no_Ns(nt_seq):
    """ Returns True if a sequence does not have any N's """

    if 'N' not in nt_seq:
        return True
    else:
        return False


def _check_vdj_len(nt_seq):
    """ Verifies the length of a VDJ nucleotide sequence mod 3 equals 1
        (The constant region continues the reading frame)
    """

    if len(nt_seq) % 3 == 1:
        return True
    else:
        return False


def _check_cdr12_lens(row, species):
    """ Verifies the length of the CDR1 and CDR2 amino acid sequences """

    # lookup expected lengths based on V call for the given species
    expected_lengths = cdr_12_aa_lengths[species][row.v_call.split('-')[0]]

    actual_lengths = (len(row.cdr1aa), len(row.cdr2aa))
    if actual_lengths in expected_lengths:
        return True
    else:
        return False


def df_vdj_qc(frame, species, verbose=False):
    """ Convenience wrapper that verifies VDJ sequences in a dataframe
        Operations along rows:
            - Fill nucleotides missing from the beginning of the V gene or end
            of the J gene and remove gap characters arising from a deletion
            - Verify the VDJ sequence is the correct length
            - Translate VDJ sequence
            - Translate CDR1 and CDR2 sequences
            - Verify the CDR1 and CDR2 sequences are the expected lengths
            - Verify the amino acid CDR3 sequence is in the translated sequence
            - Verify there are no stop codons

        Args:
            frame (pd.DataFrame): parsed change-o output
                (see parsing.load_changeo_igblast_makedb)
            species (str): 'human' or 'mouse'
            verbose (bool): print details of sequence quality control
        Returns:
            (pd.DataFrame): copy of input dataframe with additional columns:
                - (object) nt_vdj
                - (object) aa_vdj
                - (object) cdr1aa
                - (object) cdr2aa
                - (bool) ok_vdj_len
                - (bool) ok_stop_codon
                - (bool) ok_Ns
                - (bool) ok_cdr12_lens
                - (bool) ok_cdr3_in_trans
                - (bool) ok_all - AND of all other `ok_` columns
    """

    df = frame.copy()

    # load V and J gene reference sequences
    vj_refs = {}
    for gene in ['v', 'j']:
        db_path = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                '..',
                                                'data/germlines'))
        ref_fa = '{}/imgt_{}_ig_{}.fasta'.format(db_path, species, gene)
        vj_refs.update(parsing.load_gene_fasta(ref_fa))

    # append missing 'G' to two human J gene reference seqs
    if species == 'human':
        for gene in ['IGHJ6*02', 'IGHJ6*03']:
            if not vj_refs[gene].endswith('G'):
                vj_refs[gene] += 'G'

    df['nt_vdj'] = df.apply(fill_clean_sequence,
                            axis=1,
                            vj_ref_dict=vj_refs,
                            verbose=verbose)

    df['ok_vdj_len'] = df.nt_vdj.apply(_check_vdj_len)

    df['aa_vdj'] = df.nt_vdj.apply(translate)

    # check for missing CDR sequences
    for col in ['cdr1_imgt', 'cdr2_imgt', 'cdr3aa']:
        if pd.isnull(df[col]).any():
            print('\nWarning: missing {} sequence(s).'
                  ' Replacing with empty string'.format(col))
            df[col] = df[col].fillna('')

    df['cdr1aa'] = df.cdr1_imgt.apply(lambda x: translate(x.replace('.', '')))

    df['cdr2aa'] = df.cdr2_imgt.apply(lambda x: translate(x.replace('.', '')))

    df['ok_cdr12_lens'] = df.apply(_check_cdr12_lens, axis=1, species=species)

    df['ok_stop_codon'] = df.aa_vdj.apply(_no_stop_codon)

    df['ok_Ns'] = df.nt_vdj.apply(_no_Ns)

    df['ok_cdr3_in_trans'] = df.apply(lambda x: x.cdr3aa in x.aa_vdj, axis=1)

    df['ok_all'] = (df.ok_vdj_len &
                    df.ok_cdr12_lens &
                    df.ok_stop_codon &
                    df.ok_Ns &
                    df.ok_cdr3_in_trans)

    if verbose:
        if df.ok_all.all():
            print('\nAll sequences pass quality control')
        else:
            if not df.ok_vdj_len.all():
                print('\nWarning: not all VDJ nucleotide sequences have the'
                      ' correct length. See the `ok_vdj_len` dataframe'
                      ' column.')
            if not df.ok_cdr12_lens.all():
                print('\nWarning: not all CDR1 and/or CDR2 amino acid'
                      ' sequences are the correct length. See the'
                      ' `ok_cdr12_lens` dataframe column.')
            if not df.ok_stop_codon.all():
                print('\nWarning: stop codon(s) found in one or more'
                      ' sequences. See the `ok_stop_codon` dataframe'
                      ' column.')
            if not df.ok_cdr3_in_trans.all():
                print('\nWarning: the CDR3 sequence was not found in the VDJ'
                      ' amino acid sequence. See the `ok_cdr3_in_trans`'
                      ' dataframe column.')

    return df
