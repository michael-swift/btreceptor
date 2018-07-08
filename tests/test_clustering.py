from btreceptor import clustering
import pandas as pd
import numpy as np
from numpy.testing import assert_array_equal
from pandas.testing import assert_frame_equal


class TestClustering(object):

    def setup_method(self):

        self.df = pd.DataFrame([
            ['cell1', 'ARDGAAAAAL', 10, 'IGHV3-23', 'IGHJ6'],
            ['cell2', 'ARDGAAAAAK', 10, 'IGHV3-23', 'IGHJ6'],
            ['cell3', 'VMDGAAAARM', 10, 'IGHV3-23', 'IGHJ6']
        ],
                               columns=['cell', 'cdr3aa', 'cdr3aa_len',
                                        'v_call_no_allele', 'j_call_no_allele']
                               )

    def test_pw_edit(self):

        assert_array_equal(clustering.df_pw_edit(self.df),
                           np.array([[0, 0.1, 0.4],
                                     [0.1, 0, 0.4],
                                     [0.4, 0.4, 0]]))

    def test_lins_from_subset(self):

        groups = clustering.df_lineages_from_subset(self.df, 0.85)

        # cell1 and cell2 are 90% similar and belong in the same lineage
        # cell3 is in its own lineage
        assert_array_equal(np.array([0, 0, 1]), groups)

    def test_create_lineages(self):

        df_with_lins = clustering.df_add_lineages(self.df, 0.85)

        expected = self.df.copy()
        expected['lineage'] = [0, 0, 1]

        assert_frame_equal(df_with_lins, expected)
