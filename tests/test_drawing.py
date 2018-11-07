from btreceptor import drawing
import pandas as pd
import pytest
import os
import warnings


class TestDrawing(object):

    def setup_method(self):

        self.df = pd.DataFrame([
            ['cell1', 'blue', 'square', 1, 0.5, 1],
            ['cell2', 'blue', 'square', 2, 0.5, 1],
            ['cell3', 'red', 'triangle', 2, 0.5, 2],
            ['cell4', 'red', 'triangle', 2, 0.5, 2]
        ], columns=['cell', 'node_color', 'node_shape',
                    'node_size', 'node_stroke', 'lineage']
                               )

        # 2 lineages, each with two cells and one anchor 'germline' node
        self.expected_dict = {0: {'ancestor': None,
                                  'color': 'k',
                                  'shape': 'circle',
                                  'size': 1,
                                  'stroke': 0},
                              1: {'ancestor': 0,
                                  'color': 'blue',
                                  'shape': 'square',
                                  'size': 1,
                                  'stroke': 0.5},
                              2: {'ancestor': 0,
                                  'color': 'blue',
                                  'shape': 'square',
                                  'size': 2,
                                  'stroke': 0.5},
                              3: {'ancestor': None,
                                  'color': 'k',
                                  'shape': 'circle',
                                  'size': 1,
                                  'stroke': 0},
                              4: {'ancestor': 3,
                                  'color': 'red',
                                  'shape': 'triangle',
                                  'size': 2,
                                  'stroke': 0.5},
                              5: {'ancestor': 3,
                                  'color': 'red',
                                  'shape': 'triangle',
                                  'size': 2,
                                  'stroke': 0.5}}

    def test_generate_node_dict(self):

        output = drawing.df_generate_node_dict(self.df)

        assert self.expected_dict == output

    def test_drawing(self, tmpdir):

        pytest.importorskip("graph_tool")

        # missing and/or outdated modules / dependencies
        warnings.simplefilter('ignore', DeprecationWarning)
        warnings.simplefilter('ignore', RuntimeWarning)

        tmpfile = str(tmpdir.join('img.png'))

        drawing.draw_gviz(self.expected_dict, output=tmpfile)

        assert os.path.exists(tmpfile)
