from btreceptor import drawing
import pandas as pd
import pytest


def test_generate_node_dict():

    df = pd.DataFrame([
        ['cell1', 'blue', 'square', 1, 0.5, 1],
        ['cell2', 'blue', 'square', 2, 0.5, 1],
        ['cell3', 'red', 'triangle', 2, 0.5, 2],
        ['cell4', 'red', 'tranglee', 2, 0.5, 2]
    ], columns=['cell', 'node_color', 'node_shape',
                'node_size', 'node_stroke', 'lineage']
                      )

    # 2 lineages, each with two cells and one anchor 'germline' node
    expected_dict = {0: {'ancestor': None,
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
                         'shape': 'tranglee',
                         'size': 2,
                         'stroke': 0.5}}

    output = drawing.df_generate_node_dict(df)

    assert expected_dict == output
