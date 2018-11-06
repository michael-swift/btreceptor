from btreceptor import clustering
import pandas as pd
import numpy as np
from numpy.testing import assert_array_equal
from pandas.testing import assert_frame_equal


class TestClustering(object):

    def setup_method(self):

        cols = ['cell', 'cdr3aa', 'cdr3aa_len', 'v_call_no_allele',
                'j_call_no_allele']

        self.df_group_1 = pd.DataFrame([
            ['cell1', 'ARDGAAAAAL', 10, 'IGHV3-23', 'IGHJ6'],
            ['cell2', 'ARDGAAAAAK', 10, 'IGHV3-23', 'IGHJ6'],
            ['cell3', 'VMDGAAAARM', 10, 'IGHV3-23', 'IGHJ6']
        ],
                                            columns=cols)

        df_other_seqs = pd.DataFrame([
            ['cell4', 'AGLK', 4, 'IGHV4-34', 'IGHJ4'],
            ['cell5', 'ARGMMMMK', 8, 'IGHV4-34', 'IGHJ4']
        ],
                                     columns=cols)

        # ignore_index=False for testing index with duplicate values
        self.df_mult_groups = self.df_group_1.append(df_other_seqs,
                                                     ignore_index=False)

    def test_pw_edit(self):
        """ First three cells belong to same V/J/CDR3 length group """

        assert_array_equal(clustering.df_pw_edit(self.df_group_1),
                           np.array([[0, 0.1, 0.4],
                                     [0.1, 0, 0.4],
                                     [0.4, 0.4, 0]]))

    def test_lins_from_subset(self):
        """ First three cells belong to same V/J/CDR3 length group """

        groups = clustering.df_lineages_from_subset(self.df_group_1, 0.85)

        # cell1 and cell2 are 90% similar and belong in the same lineage
        # cell3 is in its own lineage
        assert_array_equal(np.array([0, 0, 1]), groups)

    def test_create_lineages(self):

        df_with_lins = clustering.df_add_lineages(self.df_mult_groups, 0.85)

        expected = self.df_mult_groups.reset_index(drop=True)
        expected['lineage'] = [0, 0, 1, 2, 3]

        assert_frame_equal(df_with_lins, expected)
