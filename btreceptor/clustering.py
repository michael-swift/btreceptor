from __future__ import division
import pandas as pd
import numpy as np
import Levenshtein
from scipy.spatial.distance import squareform
from scipy.sparse.csgraph import connected_components
from itertools import combinations


def df_pw_edit(frame):
    """ Returns array of pairwise edit distances in square form """

    ed = np.zeros(int((frame.shape[0]/2)*(frame.shape[0]-1)), dtype='float')
    for c, (x, y) in enumerate(combinations(frame.cdr3aa.values, 2)):
        ed[c] = Levenshtein.distance(x, y) / np.max([len(x), len(y)])
    sq = squareform(ed)

    return sq


def df_lineages_from_subset(frame, similarity_cutoff):
    """ Returns an array of lineage membership based on a CDR3 cutoff """

    edit_sq = df_pw_edit(frame)

    n_groups, labels = connected_components(edit_sq <= (1 - similarity_cutoff))

    return labels


def df_add_lineages(dfin, similarity_cutoff):
    """ Returns input dataframe with additional lineage column

    Args:
        similarity_cutoff (float): e.g. 0.8 for 80% minimum cdr3aa similarity
    """

    dfin = dfin.copy()

    lincnt = 0
    lins = []

    for (v, j, _), sub in dfin.groupby(['v_call_no_allele',
                                       'j_call_no_allele',
                                       'cdr3aa_len']):
        if sub.shape[0] > 1:
            # CDR3 distance comparisoin
            sub_lineages = df_lineages_from_subset(sub, similarity_cutoff)
            lins += zip(sub.index, sub_lineages + lincnt)
            lincnt += sub.shape[0]
        else:
            # single sequence belongs in its own lineage
            lins.append((sub.index.values[0], lincnt))
            lincnt += 1

    # adds a "lineage" column corresponding to the lineage number for that cell
    lins = pd.DataFrame(lins, columns=['index', 'lineage']).set_index('index')
    if 'lineage' in dfin.columns:
        dfin = dfin.drop('lineage', axis=1).join(lins)
    else:
        dfin = dfin.join(lins)

    return dfin
