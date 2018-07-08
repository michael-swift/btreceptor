from btreceptor import sequences


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
    seq_id = 'test_with_gap'
    ref_v = 'ATATATGCGCGCATGCATGTCGCATAGCGCGACGTCTA'
    ref_j = 'ATTTCGCTGCAACGGATCTGAACCGGCCTTCACACATT'
    fake_cdr3 = 'ATGCAGATGGGAACCGGG'

    # missing 1 nucleotide in beginning and end, has gap in V
    seq = ref_v[1:10] + '---' + ref_v[10:-5] + fake_cdr3 + ref_j[5:-1]

    expected = ref_v[:-5] + fake_cdr3 + ref_j[5:]
    outcome = sequences.fill_clean_sequence(seq_id, ref_v, ref_j, seq)

    assert outcome == expected
