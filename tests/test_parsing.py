from btreceptor import parsing
import os


data_dir = os.path.realpath(os.path.join(os.path.realpath(__file__),
                                         '../../data'))


def test_load_gene_fasta():

    j_genes = parsing.load_gene_fasta('{}/germlines/'
                                      'imgt_human_ig_j.fasta'.format(data_dir))

    j6 = 'ATTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG'

    assert j_genes['IGHJ6*01'] == j6
