from __future__ import division
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def _detect_dataframe_cell_type(df):
    """ Detects whether cell type being analyzed is B or T cells
    """

    if df.v_call.str.contains('IG\wV').all():
        return 'B cell'
    elif df.v_call.str.contains('TR\wV').all():
        return 'T cell'
    else:
        raise ValueError('Mixed / unknown cell types found in dataframe: '
                         '{}'.format(df.v_call.str.slice(0, 2).unique()))


def rename_changeo_columns(df, revert=False):
    """ Simplifies dataframe column names.

    Args:
        df (pd.DataFrame)
        reverse (bool): whether to revert column names to original

    Returns:
        Dataframe
    """

    df = df.copy()

    rename_dict = {'sequence_id': 'seqid',
                   'cdr3_igblast_nt': 'cdr3nt',
                   'cdr3_igblast_aa': 'cdr3aa'}
    rename_dict.update({'{}_seq_length'.format(x):
                        '{}_len'.format(x) for x in ['v', 'd', 'j']})
    rename_dict.update({'{}_seq_start'.format(x):
                        '{}_start'.format(x) for x in ['v', 'd', 'j']})

    if not revert:
        # rename columns for convenience
        df.rename(columns=lambda x: x.lower(), inplace=True)
        df.rename(columns=rename_dict, inplace=True)
    else:
        # revert to original columns
        df.rename(columns={v: k for k, v in rename_dict}, inplace=True)
        df.rename(columns=lambda x: x.upper(), inplace=True)

    if 'cdr3aa' not in df.columns:
        print('\nWarning: change-o MakeDB.py was not run with "--cdr3".'
              ' Translating the amino acid CDR3 from IMGT nucleotide'
              ' sequence instead.')
        df['cdr3aa'] = df.cdr3_imgt.apply(
            lambda x: x if pd.isnull(x) else str(Seq(x).translate()))

    return df


def load_changeo_igblast_makedb(infile):
    """ Loads tab separated Change-O output into pandas dataframe. Adds
        additional convenience columns e.g. for clustering
    """

    ig = pd.read_csv(infile, sep='\t')

    ig = rename_changeo_columns(ig)

    cell_type = _detect_dataframe_cell_type(ig)
    if cell_type == 'B cell':
        ig['heavy'] = ig.v_call.str.startswith('IGH')
    elif cell_type == 'T cell':
        ig['heavy'] = (ig.v_call.str.startswith('TRB') |
                       ig.v_call.str.startswith('TRD'))

    # update / add columns
    ig['j_end'] = ig.j_start + ig.j_len
    ig['j_start_vdj'] = ig.j_start - ig.v_start
    ig['v_call_no_allele'] = ig.v_call.str.split('*', expand=True)[0]
    ig['j_call_no_allele'] = ig.j_call.str.split('*', expand=True)[0]
    ig['cdr3aa_len'] = ig.cdr3aa.str.len()

    # cast columns as bool if not already
    bool_cols = ['functional', 'indels', 'stop', 'in_frame']
    for col in bool_cols:
        if not pd.api.types.is_bool_dtype(ig[col]):
            ig[col] = ig[col] == 'T'

    return ig


def load_gene_fasta(infile):
    """ Reads fasta into dictionary
    """

    gene_dict = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
    # convert SeqRecord sequences to uppercase string
    gene_dict = {k: str(v.seq).upper() for k, v in gene_dict.items()}

    return gene_dict
