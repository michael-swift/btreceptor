from btreceptor import sequences


def test_fill_missing_nt():
    ref_seq = 'ATATGGGGCCTCTAC'
    test_seq =  'TATGGGGCCCTAC'

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
