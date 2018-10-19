from btreceptor import sequences
import pandas as pd


def test_fill_missing_nt():
    ref_seq = 'ATATGGGGCCTCTAC'
    test_seq =  'TATGGGGCCCT'

    filled = sequences._fill_missing_nt(ref_seq, test_seq)

    assert filled == ref_seq


def test_v_fill():
    ref_seq = 'CGATTACACTGATCAGGGGTACAC'
    seq =       'ATTACACTGATCACCCGTACAC'

    filled = sequences._vfill(ref_seq, seq)

    assert filled == 'CG' + seq


def test_j_fill():
    ref_seq = 'CCATACACGATTACACTGACGG'
    seq =     'CCATAGACGATTGCACTGA'

    filled = sequences._jfill(ref_seq, seq)

    assert filled == seq + 'CGG'


def test_fill_clean():

    refs = {'IGHVX-Y': 'ATATATGCGCGCATGCATGTCGCATAGCGCGACGTCTA',
            'IGHJZ': 'ATTTCGCTGCAACGGATCTGAACCGGCCTTCACACATT'}

    # create sequence missing 1 nucleotide in beginning and end, has gap in V
    fake_cdr3 = 'ATGCAGATGGGAACCGGG'
    seq = (refs['IGHVX-Y'][1:10] +
           '---' +
           refs['IGHVX-Y'][10:-5] +
           fake_cdr3 +
           refs['IGHJZ'][5:-1])

    series = pd.Series({'seqid': 'test_with_gap',
                        'sequence_vdj': seq,
                        'v_call': 'IGHVX-Y',
                        'j_call': 'IGHJZ'
                        })

    expected = refs['IGHVX-Y'][:-5] + fake_cdr3 + refs['IGHJZ'][5:]
    outcome = sequences.fill_clean_sequence(series, refs, verbose=False)

    assert outcome == expected


def test_no_stop_codon():

    assert sequences._no_stop_codon('ASGY')


def test_stop_codon():

    assert not sequences._no_stop_codon('ASG*Y')


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
