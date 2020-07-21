from btreceptor import parsing
import os


data_dir = os.path.realpath(os.path.join(os.path.realpath(__file__),
                                         '../../data'))


def test_load_gene_fasta():

    j_genes = parsing.load_gene_fasta('{}/germlines/'
                                      'imgt_human_ig_j.fasta'.format(data_dir))

    j6 = 'ATTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG'

    assert j_genes['IGHJ6*01'] == j6


def test_parse_changeo_db_pass():
    """ Simple test for parsing an MakeDB.py output test file from
    https://bitbucket.org/kleinstein/changeo/src/3a70c5fda2a27e2ea4613eb1ce9e307b8bc06a4f/tests/data/imgt_ig_db-pass.tsv
    """

    test_data_dir = os.path.dirname(os.path.realpath(__file__))

    df = parsing.load_changeo_igblast_makedb(
        '{}/imgt_ig_db-pass.tsv'.format(test_data_dir)
    )

    assert "cdr3aa" in df.columns
