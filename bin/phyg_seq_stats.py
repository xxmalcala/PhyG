#!/usr/bin/env python


"""
Generates summary statistics and simple figures.

pip install codon-bias

Dependencies: BioPython, Matplotlib, Numpy, Seaborn, Pandas
"""
import pandas as pd
import numpy as np

from collections import defaultdict

import codonbias as cb
from Bio import SeqIO


def calc_gc(
        seq: str) -> float:
    return 100 * (seq.lower().count('g') + seq.lower().count('c')) / len(seq)


def gc_stats_per_seq(
        seq: str,
        gen_code: str = '1') -> None:

    degen_4fold_cdns = ['TC','CT','CC','CG','AC','GT','GC','GG']

    if gen_code in ['6','29']:
        degen_4fold_cdns.append('TA')

    global_gc = calc_gc(seq)
    gc1 = calc_gc("".join([seq[n] for n in range(0, len(seq), 3) if len(seq[n:n+3]) == 3]))
    gc2 = calc_gc("".join([seq[n+1] for n in range(0, len(seq), 3) if len(seq[n:n+3]) == 3]))
    gc3 = calc_gc("".join([seq[n+2] for n in range(0, len(seq), 3) if len(seq[n:n+3]) == 3]))
    gc12 = calc_gc("".join([seq[n:n+2] for n in range(0, len(seq), 3) if len(seq[n:n+3]) == 3]))

    gc3_4fld_dgn = calc_gc("".join([seq[n+2] for n in range(0, len(seq), 3) if len(seq[n:n+3]) == 3 and seq[n:n+2] in degen_4fold_cdns]))

    return global_gc, gc1, gc2, gc3, gc12, gc3_4fld_dgn


def wright_enc_exp_calc(
        gc3_4fld_dgn: float,
        null_vals: bool = False) -> float:
    if not null_vals:
        gc3_4fld = gc3_4fld_dgn/100
    else:
        gc3_4fld = gc3_4fld_dgn
    exp_enc = 2 + gc3_4fld + (29 / ((gc3_4fld**2) + (1 - gc3_4fld)**2))
    return f'{exp_enc:.2f}'


def null_ENc_GC3() -> pd.DataFrame:
    null = [float(wright_enc_exp_calc(n, True)) for n in np.arange(0,.51,0.01)]
    null += null[:-1][::-1]
    null_df = pd.DataFrame.from_dict({i:j for i, j in zip([n for n in range(0, 101)],null)}, orient = 'index')
    null_df.index.names = ['GC3-4-Fold-Degen']
    null_df = null_df.rename(columns = {0:'ENc'})
    return null_df


def get_seq_comp_enc(
        fasta_file: str,
        out_tsv: str,
        gen_code: str = '1',
        delim: str = '|'):

    enc_calc = cb.scores.EffectiveNumberOfCodons(
                            bg_correction = True,
                            genetic_code = int(gen_code))

    seq_info = defaultdict(dict)
    for i in SeqIO.parse(fasta_file,'fasta'):
        global_gc, gc1, gc2, gc3, gc12, gc3_4fld_dgn = gc_stats_per_seq(f'{i.seq}', gen_code)
        seq_info[i.id]['Length'] = len(i.seq)
        if 'Cov_' in i.id:
            seq_info[i.id]['Coverage'] = float(i.id.partition("_Cov_")[-1].partition(delim)[0])

        seq_info[i.id]['GC-Content'] = f'{global_gc:.2f}'
        seq_info[i.id]['GC1'] = f'{gc1:.2f}'
        seq_info[i.id]['GC2'] = f'{gc2:.2f}'
        seq_info[i.id]['GC3'] = f'{gc3:.2f}'
        seq_info[i.id]['GC12'] = f'{gc12:.2f}'
        seq_info[i.id]['GC3-4-Fold-Degen'] = f'{gc3_4fld_dgn:.2f}'
        seq_info[i.id]['Expected-ENc'] = wright_enc_exp_calc(gc3_4fld_dgn)
        seq_info[i.id]['ENc'] = f'{enc_calc.get_score(str(i.seq)):.2f}'

    df = pd.DataFrame.from_dict(seq_info, orient = 'index')
    df.index.names = ['Sequence']

    df.to_csv(out_tsv, sep = '\t')
