from btreceptor import sequences, parsing
import pandas as pd
import numpy as np
from numpy.testing import assert_array_equal
import os
import warnings
from Bio import BiopythonWarning
import pytest


def test_fill_missing_nt():
    ref_seq_aln = 'ATATGGGGCCTCTAC'
    seq_aligned = '--TGGGGCCCTCTAC'

    filled = sequences._fill_missing_nt(ref_seq_aln, seq_aligned)

    assert filled == ref_seq_aln[:2] + seq_aligned[2:]


def test_no_fill_necessary():
    ref_seq_aln = 'ATATGGGGCCTCTAC'
    seq_aligned = 'ATATGGGGCCTCTAC'

    filled = sequences._fill_missing_nt(ref_seq_aln, seq_aligned)

    assert filled == ref_seq_aln


def test_no_trim_necessary():
    ref_seq_aln = 'ATATGGGGCCTCTAC'
    seq_aligned = 'ATATGGGGCCTCTAC'

    filled = sequences._trim_extra_nt(ref_seq_aln, seq_aligned)

    assert filled == ref_seq_aln


def test_replace_Ns():
    ref_seq_aln = 'ATATGGGGCCTCTAC'
    seq_aligned = 'NTATGGGGNCTCTNC'

    replaced = sequences._replace_Ns_with_ref(ref_seq_aln, seq_aligned)

    assert replaced == ref_seq_aln


def test_no_replace_Ns():
    ref_seq_aln = 'ATATGGGGCCTCTAC'
    seq_aligned = 'ATATGGGGCCTCTAC'

    replaced = sequences._replace_Ns_with_ref(ref_seq_aln, seq_aligned)

    assert replaced == ref_seq_aln


class TestVgeneCorrection(object):

    def setup_method(self):

        self.ref_seq = 'CGATTACACTGATCAGGGGTACAC'

    def test_v_fill(self):
        # missing CG at the beginning
        seq = 'ATTACACTGATCACCCGTACAC'

        filled = sequences._vfix(self.ref_seq, seq)

        assert filled == 'CG' + seq

    def test_v_trim(self):
        # extra ATGC at the beginning
        seq = 'ATGCCGATTACACTGATCACCCGTACAC'

        filled = sequences._vfix(self.ref_seq, seq)

        assert filled == seq[4:]

    def test_v_retain_deletion(self):

        seq = self.ref_seq[:6] + self.ref_seq[9:]

        filled = sequences._vfix(self.ref_seq, seq)

        assert filled == seq

    def test_v_retain_insertion(self):

        seq = self.ref_seq[:6] + 'CAT' + self.ref_seq[6:]

        filled = sequences._vfix(self.ref_seq, seq)

        assert filled == seq


class TestJgeneCorrection(object):

    def setup_method(self):

        self.ref_seq = 'CCATACACGATTACACTGACGG'

    def test_j_fill(self):
        # missing final 3 nts
        seq = 'CCATAGACGATTGCACTGA'

        filled = sequences._jfix(self.ref_seq, seq)

        assert filled == seq + 'CGG'

    def test_j_trim(self):
        # extra 3 nts
        seq = 'CCATAGACGATTGCACTGACGGACG'

        filled = sequences._jfix(self.ref_seq, seq)

        assert filled == seq[:-3]

    def test_j_retain_deletion(self):

        seq = self.ref_seq[:3] + self.ref_seq[6:]

        filled = sequences._jfix(self.ref_seq, seq)

        assert filled == seq

    def test_j_retain_insertion(self):

        seq = self.ref_seq[:3] + 'CAT' + self.ref_seq[3:]

        filled = sequences._jfix(self.ref_seq, seq)

        assert filled == seq


def test_fill_clean():

    refs = {'IGHVX-Y': 'ATATATGCGCGCATGCATGTCGCATAGCGCGACGTCTA',
            'IGHJZ': 'ATTTCGCTGCAACGGATCTGAACCGGCCTTCACACATT'}

    # create sequence with 1 nucleotide missing in beginning, 2 extra
    # nucleotides at the end, a V insertion, J deletion, gap, and N
    fake_cdr3 = 'ATGCAGATGGGAACCGGG'
    seq = (refs['IGHVX-Y'][1:5] +
           'N' +
           refs['IGHVX-Y'][6:10] +
           '---' +
           refs['IGHVX-Y'][10:20] +
           'CAT' +
           refs['IGHVX-Y'][20:-5] +
           fake_cdr3 +
           refs['IGHJZ'][5:20] +
           refs['IGHJZ'][23:] + 'AT'
           )

    series = pd.Series({'seqid': 'test_with_gap',
                        'sequence_vdj': seq,
                        'v_call': 'IGHVX-Y',
                        'j_call': 'IGHJZ'
                        })

    expected = (refs['IGHVX-Y'][:20] + 'CAT' + refs['IGHVX-Y'][20:-5] +
                fake_cdr3 + refs['IGHJZ'][5:20] + refs['IGHJZ'][23:])
    outcome = sequences._fill_clean_sequence(series, refs, verbose=False)

    assert outcome == expected


def test_no_need_fill_clean():

    refs = {'IGHVX-Y': 'ATATATGCGCGCATGCATGTCGCATAGCGCGACGTCTA',
            'IGHJZ': 'ATTTCGCTGCAACGGATCTGAACCGGCCTTCACACATT'}

    # seq has full V / J and no deletions so no fill / clean necessary
    seq = refs['IGHVX-Y'] + refs['IGHJZ']

    series = pd.Series({'seqid': 'test_ok',
                        'sequence_vdj': seq,
                        'v_call': 'IGHVX-Y',
                        'j_call': 'IGHJZ'
                        })

    expected = seq
    outcome = sequences._fill_clean_sequence(series, refs, verbose=False)

    assert outcome == expected


def test_no_stop_codon():

    assert sequences._no_stop_codon('ASGY')


def test_stop_codon():

    assert not sequences._no_stop_codon('ASG*Y')


def test_Ns_in_nt_seq():

    assert not sequences._no_Ns('ATGNNC')


def test_no_Ns_in_nt_seq():

    assert sequences._no_Ns('ATGC')


def test_correct_vdj_len():

    assert sequences._check_vdj_len('CAGG')


def test_wrong_vdj_len():

    assert not sequences._check_vdj_len('CAG')
    assert not sequences._check_vdj_len('CAGGG')


def test_translate():

    assert sequences.translate('GAGGTGCAACTG') == 'EVQL'


def test_cdr12_correct_lens():

    species = 'human'
    series = pd.Series({'v_call': 'IGHV1',
                        'cdr1aa': 'GLPSNDNW',
                        'cdr2aa': 'ISPYDSNT'})

    assert sequences._check_cdr12_lens(series, species)


def test_cdr12_wrong_lens():

    species = 'human'
    series = pd.Series({'v_call': 'IGHV1',
                        'cdr1aa': 'GL',
                        'cdr2aa': 'IS'})

    assert not sequences._check_cdr12_lens(series, species)


class TestVDJqc(object):

    def setup_method(self):
        # load test immune repertoire data
        test_data_dir = os.path.dirname(os.path.realpath(__file__))

        self.df = parsing.load_changeo_igblast_makedb(
            '{}/imgt_ig_db-pass.tsv'.format(test_data_dir))

    def test_vdj_qc(self):

        # convert v/j_call columns to genes
        # e.g. Homsap IGLV1-40*01 F,Homsap IGLV1-40*02 F --> IGLV1-40*01
        self.df['v_call'] = self.df.v_call.apply(lambda x: x.split(' ')[1])
        self.df['j_call'] = self.df.j_call.apply(lambda x: x.split(' ')[1])

        warnings.simplefilter('ignore', BiopythonWarning)  # skip partial codon
        dfqc = sequences.df_vdj_qc(self.df, 'human', verbose=True)

        # first 3 sequences have N's, otherwise remaining 4 are OK
        ok_Ns = np.array([False]*3 + [True]*4)

        assert_array_equal(dfqc.ok_Ns.values, ok_Ns)
        assert_array_equal(dfqc.ok_all.values, ok_Ns)

    def test_vdj_unexpected_gene_calls(self):

        # don't correct v_call / j_call, expect KeyError
        with pytest.raises(KeyError):
            sequences.df_vdj_qc(self.df, 'human', verbose=True)
