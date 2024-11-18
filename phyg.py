#!/usr/bin/env python3

"""
Controls the flow for PhyG's modules. Currently does NOT resume in-progress
runs. This will be fixed. Additionally, there will be an "adding" data suite of
modules/scripts, that are yet to come [as of 2024-05-11].
"""

import logging,sys, time
from datetime import datetime, timedelta
from pathlib import Path


from bin import phyg_eval_args as eval_args


def check_args():
    args = eval_args.collect_args()

    return args


# def log_info(qparams: list, arg_vars: dict) -> None:
#     for k in qparams:
#         x = f'{k}:'.ljust(22)
#         logging.info(f'{x}{str(arg_vars[k.replace("-","_")] or "None")}')
#
# def init_log_stats(args):
#     arg_vars = vars(args)
#
#     curr_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
#
#     initdir = f'{arg_vars["msa_dir"]}'
#     outdir = f'{arg_vars["project_name"]}_PhyG'
#
#     Path(outdir).mkdir(parents = True, exist_ok = True)
#
#     outlog = f'{outdir}/tides.log'
#
#     logging.basicConfig(filename = outlog,
#         filemode = 'w+',
#         format = '%(message)s',
#         level = logging.DEBUG)
#
#     logging.info(
#         f'+---------------------------------+\n'\
#         f'|         PhyG Parameters         |\n'\
#         f'+---------------------------------+'
#         )
#
#     logging.info(f'\nRun-Start:            {curr_time}')
#     logging.info(f'Query-Dir:            {initdir.rstrip("/")}')
#     logging.info(f'Output-Dir:           {outdir}')
#
#     if arg_vars['no_filtering']:
#         qparams = ['threads']
#         log_info(qparams, arg_vars)
#         logging.info(f'Homology-Filtering:       NONE')
#
#     elif arg_vars['guidance']:
#         qparams = ['max-iter','seq-score','threads']
#         log_info(qparams, arg_vars)
#         logging.info(f'Homology-Filtering:   Guidance2')
#
#     elif arg_vars['iqr']:
#         qparams = ['threads']
#         log_info(qparams, arg_vars)
#         logging.info(f'Homology-Filtering:   P-dist inter-quartile range')
#
#     else:
#         qparams = ['threads']
#         log_info(qparams, arg_vars)
#         logging.info(f'Homology-Filtering:   Bootstrap p-dist')
#
#     return outlog


if __name__ == '__main__':

    # args = check_args()
    print('Currently re-working argument handling...')

    print('Thanks for using PhyG!')
